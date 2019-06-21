function [glb_res,IVEC,JVEC,VVEC,count] = cv_contribution(coords,nd_vol,iblank,nd_dof_map,sol_np1,sol_n,sol_nm1,time_info,pp, ...
                                                          glb_res,IVEC,JVEC,VVEC,count)
% control volume contributions

% loop over all nodes
for nd = 1:size(coords,1)

    % skip loop for fringe nodes
    if( iblank(nd) == -1 )
        continue;
    end

    % extract only coordinates for the current node
    coord = coords(nd,:);

    % extract dof ids associated with current node
    nd_dofs = nd_dof_map(nd,:)';

    % compute residual and jacobian contribution from volumetric
    % contributions
    [res_vol, jac_vol] = compute_res_jac_vol(sol_np1(nd_dofs), sol_n(nd_dofs), sol_nm1(nd_dofs), ...
                                             coord, nd_vol(nd), time_info, pp);

    % sum into global reidual and jacobian
    glb_res(nd_dofs) = glb_res(nd_dofs) + res_vol;

    % store jac_vol in CSC format
    for r = 1:length(nd_dofs)
        for c = 1:length(nd_dofs)
            count = count + 1;
            IVEC(count,1) = nd_dofs(r);
            JVEC(count,1) = nd_dofs(c);
            VVEC(count,1) = jac_vol(r,c);
        end
    end

end

end