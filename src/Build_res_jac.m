function [glb_res,glb_jac] = Build_res_jac(glb_mesh,glb_donor_map,glb_iblank, ....
                                           glb_sol_np1,glb_sol_n,glb_sol_nm1,glb_nd_dof_map, ...
                                           ov_info,time_info,pp)
% assemble residual and jacobian correpesonding to both grids

%% pre-processing

num_grids    = ov_info('num grids'); % nmuber of grids
ndof_nd      = pp('dof per node'); % number of dofs per node
jac_srf_size = (2*ndof_nd)^2; % size of jacobain contribution from each scs

num_dofs = zeros(num_grids,1);
tot_dofs = size(glb_sol_np1,1); % total dof counter
vec_size = 0; % vector size counter
count    = 0; % counter to keep track of IVEC, JVEC, VVEC, indexing

% loop over number of grids
for ig = 1:num_grids
     
    coords = glb_mesh{ig}{2,1}; % extract coordinates
    
    edge_topo_x = glb_mesh{ig}{2,2}; % extract edge topologies
    
    edge_topo_y = glb_mesh{ig}{2,4}; % extract edge areas
    
    num_dofs(ig) = ndof_nd*size(coords,1); % number of dofs on this mesh
    
    donor_map = glb_donor_map{ig}; % donor map array for current grid
    
    % determine size for compressed sparse column format vector
    vec_size = vec_size ...
             + jac_srf_size*(size(edge_topo_x,1) + size(edge_topo_y,1)) ...
             + ndof_nd*(size(coords,1) - size(donor_map,1));
    
    % account for fringe contributions if this is a decoupled solve
    if (ov_info('solve type') == "coupled")
        vec_size = vec_size + ndof_nd*(sum(cellfun('length',donor_map(:,2))) + size(donor_map,1));
    elseif (ov_info('solve type') == "decoupled iterative")
        vec_size = vec_size + ndof_nd*size(donor_map,1);
    end
    
end

% initialize global residual array
glb_res = zeros(tot_dofs,1);

% compressed sparse column format vector components
IVEC = zeros(vec_size,1);
JVEC = zeros(vec_size,1);
VVEC = zeros(vec_size,1);

%% loop over all grids and assemble residual and jacobian
for ig = 1:num_grids
    
    % extract relevant data for this mesh grid
    coords = glb_mesh{ig}{2,1}; % extract coordinates 

    edge_topo_x = glb_mesh{ig}{2,2}; % x-edge topology
    edge_area_x = glb_mesh{ig}{2,3}; % x-edge area

    edge_topo_y = glb_mesh{ig}{2,4}; % y-edge topology
    edge_area_y = glb_mesh{ig}{2,5}; % y-edge area

    nd_vol = glb_mesh{ig}{2,6}; % extract nodal volumes
    
    nd_dof_map = glb_nd_dof_map{ig}; % extract node to dof map
    
    iblank = glb_iblank{ig}; % iblank array for current grid

    donor_map = glb_donor_map{ig}; % donor map array for fringe nodes

    % extract donor mesh info
    donor_grid = ov_info(strcat('mesh',num2str(ig),' donor'));
    
    donor_coords = glb_mesh{donor_grid}{2,1}; % extract donor mesh coordinates
    
    donor_nd_dof_map = glb_nd_dof_map{donor_grid}; % extract node to dof map
    
    % loop over all edges aligned with x for surface contributions
    [glb_res,IVEC,JVEC,VVEC,count] = internal_scs_contribution(1,edge_topo_x,edge_area_x, ...
                                                               coords,iblank,nd_dof_map,glb_sol_np1, ...
                                                               pp, ...
                                                               glb_res,IVEC,JVEC,VVEC,count); 
    
    % loop over all edges aligned with y for surface contributions
    [glb_res,IVEC,JVEC,VVEC,count] = internal_scs_contribution(2,edge_topo_y,edge_area_y, ...
                                                               coords,iblank,nd_dof_map,glb_sol_np1, ...
                                                               pp, ...
                                                               glb_res,IVEC,JVEC,VVEC,count); 
                                                   
    % loop over all nodes for volume contribution
    [glb_res,IVEC,JVEC,VVEC,count] = cv_contribution(coords,nd_vol,iblank,nd_dof_map, ...
                                                     glb_sol_np1,glb_sol_n,glb_sol_nm1, ...
                                                     time_info,pp, ...
                                                     glb_res,IVEC,JVEC,VVEC,count);

     % skip contribution from fringe points for fully decoupled solve
     if (ov_info('solve type') == "decoupled" && ov_info('fringe update') == "direct")
         continue
     end
     
    % loop over all fringe nodes for fringe interpolation contribution
    [glb_res,IVEC,JVEC,VVEC,count] = fringe_contribution(ov_info,donor_map, ...
                                                         coords,donor_coords,nd_dof_map,donor_nd_dof_map, ...
                                                         glb_sol_np1, ...
                                                         glb_res,IVEC,JVEC,VVEC,count,false);
    
end

% create sparse matrix out of IVEC,JVEC,VVEC
glb_jac = sparse(IVEC,JVEC,VVEC,tot_dofs,tot_dofs);

end