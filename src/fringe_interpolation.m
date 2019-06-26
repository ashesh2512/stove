function sol_np1 = fringe_interpolation(glb_mesh,glb_donor_map,sol_np1,glb_nd_dof_map,ov_info)
% fringe interpolation

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
    
    % loop over all fringe nodes
    for ifr = 1:size(donor_map,1)

        frng_nd   = donor_map{ifr,1}; % fringe node
        donor_nds = donor_map{ifr,2}; % extract donor node ids

        frng_coords     = coords(frng_nd,:); % extract fringe node coordinates
        donor_nd_coords = donor_coords(donor_nds,:); % extract donor node coordinates

        frng_nd_dofs  = reshape(nd_dof_map(frng_nd,:)',[],1); % extract fringe node dofs
        donor_nd_dofs = reshape(donor_nd_dof_map(donor_nds,:)',[],1); % extract donor node dofs 

        % compute interpolation coefficients
        coeff = compute_frg_coeff(frng_coords,donor_nd_coords,ov_info);

        % determine fringe point solution based on coefficients computed above
        for r = 1:length(frng_nd_dofs)
            % interpolate solution at fringe point = N*don_sol
            sol_np1(frng_nd_dofs(r)) = coeff'*sol_np1(donor_nd_dofs);
        end
    end
    
end

end