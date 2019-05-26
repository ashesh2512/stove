% 2D overset solver using structured vertex-centered FVM.
% 
% Assumptions - non-dirichlet boundary conditions cannot be applied
%             - mesh domain alignes with x and y axis 
%             - each boundary is associated with a BC

%% clear everything
clc; close all; clear all;

%% user input options

% elemental properties
prblm       = 1; % problem id
ndof_per_nd = 1; % number of dofs per node
k1          = 1; % conductivity of domain 1
ep          = {[prblm,ndof_per_nd,k1]}; % element property

% background mesh
box1 = [-2,2; -2,2];
h1   = [ 0.2,  0.2]/2/2/2;

% boundary condition map - bottom, right, top, left
bc = containers.Map({   'bottom',     'right',       'top',      'left'}, ...
                    {"dirichlet", "dirichlet", "dirichlet", "dirichlet"});

% convergence parameters
N_iters    =    10; % Maximum number of Newton steps
resnrmdrop = 1e-09; % Newton convergence criteria

%% build mesh

mesh_obj1 = meshgen(box1, h1);

%% initialize solution and set boundary conditions

% build nd to dof map
coords = mesh_obj1{2,1};
nd_dof_map = zeros(size(coords,1),ndof_per_nd);
for ind = 1:size(coords,1)
    nd_dof_map(ind,:) = ((ind-1)*ndof_per_nd+1):(ind*ndof_per_nd);
end

% apply initial and boundary conditions
[glb_sol,fdof] = apply_ic_bc(ep{1},mesh_obj1,nd_dof_map);

%% perform Newton solve

% loop over newton steps
for in = 1:N_iters
    
    % evaluate residual and jacobian matrices for all grids
    [RES, JAC_MAT] = Build_res_jac(mesh_obj1, glb_sol, nd_dof_map, ep);
            
    fprintf('residual norm = %e \n', norm(RES(fdof)));
    if (norm(RES(fdof)) < resnrmdrop)
        break;
    elseif in == N_iters
        error('Maximum Newton iterations reached, and solution still not converged');
    end
    
    dsol = JAC_MAT(fdof,fdof)\RES(fdof); % perform linear solve
    
    glb_sol(fdof) = glb_sol(fdof) - dsol; % update solution
    
end

%% compute L2 error

L2_err = comp_L2_err(ep{1},glb_sol,coords);
fprintf('\n L2 error = %e', L2_err);
fprintf('\n');

%% Plot the solution

% 1D mesh along x
omega_x  = box1(1,2) - box1(1,1);
nnodes_x = floor(omega_x/h1(1)) + 1;
coord_x  = linspace(box1(1,1), box1(1,2), nnodes_x)';

% 1D mesh along y
omega_y  = box1(2,2) - box1(2,1);
nnodes_y = floor(omega_y/h1(2)) + 1;
coord_y  = linspace(box1(2,1), box1(2,2), nnodes_y)';

figure()
surf_handle = surf(coord_x,coord_y,reshape(glb_sol,[nnodes_x,nnodes_y])','edgecolor','none');
colorbar;
set(gcf,'color','w');
set(gca, 'FontSize', 18);
view(2);
