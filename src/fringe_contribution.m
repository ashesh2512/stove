function [glb_res,IVEC,JVEC,VVEC,count] = fringe_contribution(ov_info,donor_map, ...
                                                              coords,donor_coords,nd_dof_map,donor_nd_dof_map,sol_np1, ...
                                                              glb_res,IVEC,JVEC,VVEC,count)
% fringe interpolation contributions

% loop over all fringe nodes
for ifr = 1:size(donor_map,1)

    frng_nd   = donor_map{ifr,1}; % fringe node
    donor_nds = donor_map{ifr,2}; % extract donor node ids

    frng_coords     = coords(frng_nd,:); % extract fringe node coordinates
    donor_nd_coords = donor_coords(donor_nds,:); % extract donor node coordinates

    frng_nd_dofs  = reshape(nd_dof_map(frng_nd,:)',[],1); % extract fringe node dofs
    donor_nd_dofs = reshape(donor_nd_dof_map(donor_nds,:)',[],1); % extract donor node dofs 

    % compute residual and jacobian contribution from interpolation
    % contributions
    coeff = compute_frg_coeff(frng_coords,donor_nd_coords,ov_info);

    % sum into global reidual and jacobian
    for r = 1:length(frng_nd_dofs)
        % sum into residual = frng_sol - N*don_sol
        glb_res(frng_nd_dofs(r)) = glb_res(frng_nd_dofs(r)) ... 
                                 + sol_np1(frng_nd_dofs(r)) - coeff'*sol_np1(donor_nd_dofs);

        % jacobian contribution from d(frng_sol)/d(frng_sol)
        count = count + 1;
        IVEC(count,1) = frng_nd_dofs(r);
        JVEC(count,1) = frng_nd_dofs(r);
        VVEC(count,1) = 1.0;

        % jacobian contribution from d(-N*don_sol)/d(don_sol)
        for c = 1:length(donor_nd_dofs)
            count = count + 1;
            IVEC(count,1) = frng_nd_dofs(r);
            JVEC(count,1) = donor_nd_dofs(c);
            VVEC(count,1) = -coeff(c);
        end
    end

end

end