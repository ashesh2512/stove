function [glb_res,IVEC,JVEC,VVEC,count] = internal_scs_contribution(dim,edge_topo,edge_area, ...
                                                                    coords,iblank,nd_dof_map,sol,pp, ...
                                                                    glb_res,IVEC,JVEC,VVEC,count)
% sub control surface contribution perpendicular to dimension dim

% loop over all sub control surfaces perpendicular to dimension dim
for e = 1:size(edge_topo,1)

    edge_nd = edge_topo(e,:); % extract nodes belonging to the edge

    % extract dimension-specific information for the edge
    if(dim==1) % x-axis
        
        % extract only x-coordinates for the current edge
        edge_coord = coords(edge_nd,1)';
        
        % define normal to the edge in +ve x direction
        edge_nrml = [1; 0];   
        
    elseif (dim==2)  % y-axis
        
        % extract only y-coordinates for the current edge
        edge_coord = coords(edge_nd,2)';
        
        % define normal to the edge in +ve x direction
        edge_nrml = [0; 1];
        
    end

    % extract dof ids associated with nodes of current edge
    edge_dofs = reshape(nd_dof_map(edge_nd,:)',[],1);

    % extract field dofs
    field_dofs = reshape(nd_dof_map(edge_nd(iblank(edge_nd)==1),:)',[],1);

    % compute residual and jacobian contribution from surface
    % contributions
    [res_srf, jac_srf] = compute_res_jac_srf( sol(edge_dofs), edge_coord, edge_area(e), edge_nrml, pp );

    % sum into global reidual and jacobian for field dofs
    glb_res(field_dofs) = glb_res(field_dofs) + res_srf(ismember(edge_dofs,field_dofs));

    % store jac_srf in CSC format
    for r = 1:length(edge_dofs)
        for c = 1:length(edge_dofs)
            count = count + 1;
            IVEC(count,1) = edge_dofs(r);
            JVEC(count,1) = edge_dofs(c);

            % add jacobian contribution only for field points
            if( ismember(edge_dofs(r),field_dofs) ) 
                VVEC(count,1) = jac_srf(r,c);
            end
        end
    end

end

end