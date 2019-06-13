% 2D overset solver using structured vertex-centered FVM.
% 
% Assumptions - overset configuration with exactly 2 grids
%             - grids are stationary in time
%             - quad4 elements with constant spacing in each direction
%             - square/rectangular problem domain alignes with x and y axis
%             - same physics with same material properties solved on both grids
%             - non-dirichlet boundary conditions cannot be applied
%             - each boundary is associated with a BC
%             - only coupled overset solver supported currently

%% clear everything
clc; close all; clear all;
format long;

%% user input options

% problem properties
pp = containers.Map({            'prblm', 'dof per node', 'conductivity'}, ...
                     { "steady heat MMS",              1,              1});
                 
% background mesh
box1 = [-2,2; -2,2];
h1   = [ 0.2,  0.2]/2;

% boundary condition map - bottom, right, top, left
bc1 = containers.Map({   'bottom',     'right',       'top',      'left'}, ...
                     {"dirichlet", "dirichlet", "dirichlet", "dirichlet"});
% near body mesh
box2 = [-1.13625,0.86375; -1.13625,0.86375];
h2   = [             0.2,              0.2]/2;

% boundary condition map - bottom, right, top, left
bc2 = containers.Map({ 'bottom',   'right',     'top',    'left'}, ...
                     {"overset", "overset", "overset", "overset"});

% implicit hole cutting options -- apply to both grids
% mandatory firnge - number of points inward from the overset boundary to
%                    be marked as fringe
% overlap          - minimum distance between innermost fringe points on 
%                    each grid this value can vary from boundary to 
%                    boundary
ov_info = containers.Map({ 'num grids', 'mesh1 donor', 'mesh2 donor', 'mandatory frng', 'overlap', 'donor grid', 'intrp order'}, ...
                         {  2, 2, 1, 2, [3*h2(1), 3*h2(1), 3*h2(1), 3*h2(1)], "tensor", 1 });
% ov_info = containers.Map({ 'num grids', 'mesh1 donor', 'mesh2 donor', 'mandatory frng', 'overlap', 'donor grid', ...
%                            'intrp radius', 'intrp type', 'intrp shape', 'shape param', 'poly order'}, ...
%                          { 2, 2, 1, 2, [3*h2(1), 3*h2(1), 3*h2(1), 3*h2(1)], "radial", ...
%                            2.5*max(h2), "rbf", "gaussian", 1.0, 1 });

% convergence parameters
N_iters    =    10; % Maximum number of Newton steps
resnrmdrop = 1e-09; % Newton convergence criteria

% debug/display flags
plot_mesh_flag     = false;
plot_hole_cut_flag = true;
print_frng_dist    = true;
plot_sol_flag      = true;

%% build meshes

mesh_obj1 = meshgen(box1, h1);
mesh_obj2 = meshgen(box2, h2);

if (plot_mesh_flag)
    plot_mesh(mesh_obj1, mesh_obj2);
end

%% perform hole cutting for each mesh

% cut hole based on mandatory fringe boundaries for mesh 1
[nb_iblank1, bg_iblank1] = cut_hole(mesh_obj1, mesh_obj2, bc1, ov_info);

% cut hole based on mandatory fringe boundaries for mesh 2
[nb_iblank2, bg_iblank2] = cut_hole(mesh_obj2, mesh_obj1, bc2, ov_info);

% merge iblanks
iblank1 = nb_iblank1; iblank1(bg_iblank2==-1) = -1;
iblank2 = bg_iblank1; iblank2(nb_iblank2==-1) = -1;

% determine nodal donor maps for each mesh
[donor_map1, don_nds2] = map_donors(mesh_obj1, mesh_obj2, iblank1, iblank2, ov_info);
[donor_map2, don_nds1] = map_donors(mesh_obj2, mesh_obj1, iblank2, iblank1, ov_info);

if (plot_hole_cut_flag)
    plot_hole_cut(mesh_obj1, mesh_obj2, iblank1, iblank2, don_nds1, don_nds2);
end

if (print_frng_dist)
    compute_frng_gap(mesh_obj1,mesh_obj2,nb_iblank1,bg_iblank1,nb_iblank2,bg_iblank2);
end

%% initialize solution and set boundary conditions

ndof_per_nd = pp('dof per node'); % number of dofs per node

% build nd to dof map for mesh 1
coords1 = mesh_obj1{2,1};
nd_dof_map1 = zeros(size(coords1,1),ndof_per_nd);
for ind = 1:size(coords1,1)
    nd_dof_map1(ind,:) = ((ind-1)*ndof_per_nd+1):(ind*ndof_per_nd);
end

% build nd to dof map for mesh 2 accounting for dofs in mesh 1
coords2 = mesh_obj2{2,1};
nd_dof_map2 = zeros(size(coords2,1),ndof_per_nd);
for ind = 1:size(coords2,1)
    nd_dof_map2(ind,:) = max(nd_dof_map1(:)) + (((ind-1)*ndof_per_nd+1):(ind*ndof_per_nd));
end

% apply boundary conditions to mesh 1
[sol1,fdof1] = apply_bc(pp,mesh_obj1,bc1,nd_dof_map1);

% apply boundary conditions to mesh 2
[sol2,fdof2] = apply_bc(pp,mesh_obj2,bc2,nd_dof_map2);

% asemble global solution and free dof array
glb_sol  = [ sol1; sol2];
glb_fdof = [fdof1;fdof2]; 

%% perform coupled Newton solve 

% loop over newton steps
for in = 1:N_iters
    
    % evaluate residual and jacobian matrices for all grids
    [RES, JAC_MAT] = Build_res_jac_coupled({mesh_obj1,mesh_obj2}, {donor_map1,donor_map2}, {iblank1,iblank2},...
                                            glb_sol, {nd_dof_map1,nd_dof_map2}, ov_info, pp);
            
    fprintf('residual norm = %e \n', norm(RES(glb_fdof)));
    if (norm(RES(glb_fdof)) < resnrmdrop)
        break;
    elseif in == N_iters
        error('Maximum Newton iterations reached, and solution still not converged');
    end
    
    dsol = JAC_MAT(glb_fdof,glb_fdof)\RES(glb_fdof); % perform linear solve
    
    glb_sol(glb_fdof) = glb_sol(glb_fdof) - dsol; % update solution
    
end

%% plot solution

if (plot_sol_flag)
    plot_sol(mesh_obj1, mesh_obj2, glb_sol, nd_dof_map1, nd_dof_map2);
end

%% compute L2 error

comp_L2_err(pp,glb_sol,mesh_obj1,mesh_obj2);
