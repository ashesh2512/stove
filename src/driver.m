function L2_err = driver(inp_container)
% 2D overset solver using structured vertex-centered FVM.
% 
% Assumptions - overset configuration with exactly 2 grids
%             - grids are stationary in time
%             - quad4 elements with constant spacing in each direction
%             - square/rectangular problem domain alignes with x and y axis
%             - same physics with same material properties solved on both grids
%             - each boundary is associated with a BC
%             - non-dirichlet boundary conditions cannot be applied
%             - only coupled overset solver supported currently
%             - only BDF 1 and BDF 2 time solvers supported

%% extract input options

% extract problem properties
pp = inp_container('problem definition');
                 
% extract mesh 1 details
mesh1 = inp_container('mesh 1');
box1  = mesh1('dim');
h1    = mesh1('size');
bc1   = mesh1('bc');

% extract mesh 2 details
mesh2 = inp_container('mesh 2');
box2  = mesh2('dim');
h2    = mesh2('size');
bc2   = mesh2('bc');
                   
% extract overset properties
ov_info = inp_container('overset prop');

% extract time solver properties
time_info = inp_container('time solver prop');
BDF_order = time_info('BDF order');
init_time = time_info('init time');
tot_time  = time_info('total time');
dt        = time_info('time step');

% extract linear solver info
lin_sol_info = inp_container('lin solver prop');

% perform certain checks in regards to combination of solver types and
% fringe update option selected
if (ov_info('solve type') == "coupled" && isKey(ov_info,'fringe update'))
    if (ov_info('fringe update') == "iterative")
        error('Invalid fringe update for decoupled meshes using a direct solver');
    end
end
if (ov_info('solve type') == "decoupled" && lin_sol_info('type') == "direct" && isKey(ov_info,'fringe update'))
    if (ov_info('fringe update') == "iterative")
        error('Invalid fringe update for decoupled meshes using a direct solver');
    end
end

% extract debug/display flags
debug_flags        = inp_container('debug flags');
plot_mesh_flag     = debug_flags('plot mesh');
plot_hole_cut_flag = debug_flags('plot hole cut');
print_frng_dist    = debug_flags('print fringe gap');
plot_sol_flag      = debug_flags('plot sol');

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

% extract the free and constrained parts of the system for mesh 1
[fdof1,fld_dof1,dnd1] = get_dof_status(mesh_obj1,bc1,nd_dof_map1,ov_info,iblank1,lin_sol_info);

% extract the free and constrained parts of the system for mesh 2
[fdof2,fld_dof2,dnd2] = get_dof_status(mesh_obj2,bc2,nd_dof_map2,ov_info,iblank2,lin_sol_info);

% asemble global array of free dofs 
glb_fdof = [fdof1; fdof2]; 

% asemble global array of field dofs 
glb_fld_dof = [fld_dof1; fld_dof2]; 

% apply initial conditions to mesh 1
[sol1_n] = apply_ic_n(pp,mesh_obj1,init_time);

% apply initial conditions to mesh 2
[sol2_n] = apply_ic_n(pp,mesh_obj2,init_time);
    
% asemble global solution array at n
glb_sol_n = [sol1_n; sol2_n];

% initialize global solution array at n-1
glb_sol_nm1 = zeros(size(glb_sol_n,1),1);
    
if (plot_sol_flag)
    plot_sol(mesh_obj1, mesh_obj2, glb_sol_n, nd_dof_map1, nd_dof_map2, 0.0);
end

%% Time loop

time_step_count = 0; % time step counter
curr_time       = init_time; % intialize time of solution

while curr_time <= tot_time
    
    time_step_count = time_step_count+1; % update time step count
    curr_time = curr_time + dt; % update simulation time
    fprintf('\nTime = %e \n\n', curr_time);
    
    % set time stepping parameters
    if(((time_step_count == 1) && (BDF_order == 2)) || (BDF_order == 1))
        time_info('gamma_np1') = 1.0;
        time_info('gamma_n1') = -1.0;
        time_info('gamma_nm1') = 0.0; 
    elseif((time_step_count > 1) && (BDF_order == 2))
        time_info('gamma_np1') = 1.5;
        time_info('gamma_n1') = -2;
        time_info('gamma_nm1') = 0.5;
    else
        error('Invalid BDF order');
    end
    
    % apply boundary conditions to mesh 1
    sol1_np1 = apply_bc(pp,mesh_obj1,dnd1,nd_dof_map1,curr_time);

    % apply boundary conditions to mesh 2
    sol2_np1 = apply_bc(pp,mesh_obj2,dnd2,nd_dof_map2,curr_time);
    
    % asemble global solution array at n+1
    glb_sol_np1 = [sol1_np1; sol2_np1];
    
    switch ov_info('solve type')
        case "coupled" % perform strongly coupled linear solve
            [glb_sol_np1,glb_sol_n,glb_sol_nm1] = solver_coupled(mesh_obj1,mesh_obj2,donor_map1,donor_map2,iblank1,iblank2, ...
                                                                 glb_sol_np1,glb_sol_n,glb_sol_nm1,nd_dof_map1,nd_dof_map2,glb_fdof, ...
                                                                 ov_info,time_info,lin_sol_info,pp);
                                                               
        case "decoupled" % perform weakly coupled linear solve
            switch lin_sol_info('type')
                case "direct"
                    [glb_sol_np1,glb_sol_n,glb_sol_nm1] = solver_decoupled(mesh_obj1,mesh_obj2,donor_map1,donor_map2,iblank1,iblank2, ...
                                                                           glb_sol_np1,glb_sol_n,glb_sol_nm1,nd_dof_map1,nd_dof_map2,glb_fdof, ...
                                                                           ov_info,time_info,lin_sol_info,pp);
                                                                       
                case "iterative"
                    [glb_sol_np1,glb_sol_n,glb_sol_nm1] = solver_decoupled_iterative(mesh_obj1,mesh_obj2,donor_map1,donor_map2,iblank1,iblank2, ...
                                                                                     glb_sol_np1,glb_sol_n,glb_sol_nm1,nd_dof_map1,nd_dof_map2,glb_fdof,glb_fld_dof, ...
                                                                                     ov_info,time_info,lin_sol_info,pp);
            end
                                                               
        otherwise
            error('Invalid overset solve type');
    end
                   
    % plot solution
    if (plot_sol_flag)
        plot_sol(mesh_obj1, mesh_obj2, glb_sol_np1, nd_dof_map1, nd_dof_map2, curr_time);
    end 
    
end

%% compute L2 error

L2_err = comp_L2_err(pp,glb_sol_np1,mesh_obj1,mesh_obj2,curr_time);

end

