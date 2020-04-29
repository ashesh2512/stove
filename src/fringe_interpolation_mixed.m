function frg_sol = fringe_interpolation_mixed(fringe_grid,donor_grid,donor_map,frg_sol,don_sol,frg_dof_map,don_dof_map,ov_info)
% grid-specific fringe interpolation

nb_coord = fringe_grid{2,1}; % extract fringe mesh coordinates
bg_coord = donor_grid{2,1}; % extract donor mesh coordinates

% loop over all fringe nodes
for ifr = 1:size(donor_map,1)
    
    frng_pts  = donor_map{ifr,1}; % fringe node
    donor_pts = donor_map{ifr,2}; % extract donor node ids
    
    frng_coords  = nb_coord(frng_pts,:); % extract fringe node coordinates
    donor_coords = bg_coord(donor_pts,:); % extract donor node coordinates
    
    frng_dofs  = reshape(frg_dof_map(frng_pts,:)',[],1); % extract fringe node dofs
    donor_dofs = reshape(don_dof_map(donor_pts,:)',[],1); % extract donor node dofs
    
    % compute interpolation coefficients
    coeff = compute_frg_coeff(frng_coords,donor_coords,ov_info);
    
    % determine fringe point solution based on coefficients computed above
    for r = 1:length(frng_dofs)
        % interpolate solution at fringe point = N*don_sol
        frg_sol(frng_dofs(r)) = coeff'*don_sol(donor_dofs);
    end
end
    
end