function glb_res = fringe_DBC(glb_res,glb_mesh,glb_donor_map,fulldsol,glb_nd_dof_map,ov_info)
% fringe residual contribution during iterative solve

num_grids = ov_info('num grids'); % nmuber of grids

% loop over all grids and determine solution at fringe points
for ig = 1:num_grids
    
    % extract relevant data for this mesh grid
    coords = glb_mesh{ig}{2,1}; % extract coordinates 
    
    nd_dof_map = glb_nd_dof_map{ig}; % extract node to dof map
    
    donor_map = glb_donor_map{ig}; % donor map array for fringe nodes

    % extract donor mesh info
    donor_grid = ov_info(strcat('mesh',num2str(ig),' donor'));
    
    donor_coords = glb_mesh{donor_grid}{2,1}; % extract donor mesh coordinates
    
    donor_nd_dof_map = glb_nd_dof_map{donor_grid}; % extract node to dof map
   
    glb_res = fringe_contribution(ov_info,donor_map, ...
              coords,donor_coords,nd_dof_map,donor_nd_dof_map,fulldsol, ...
              glb_res,[],[],[],0,true);
    
end

end