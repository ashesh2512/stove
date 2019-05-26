function [glb_res,glb_jac] = Build_res_jac(mesh,glb_sol,nd_dof_map,ep)

% initialize global residual and jacobian arrays
tot_dofs = size(glb_sol,1);
glb_res  = zeros(tot_dofs, 1);
glb_jac  = zeros(tot_dofs, tot_dofs);

% extract coordinates 
coords = mesh{2,1};

% extract edge topologies
edge_topo_x = mesh{2,2};
edge_area_x = mesh{2,3};

% extract edge areas
edge_topo_y = mesh{2,4};
edge_area_y = mesh{2,5};

% extract nodal volumes
nd_vol = mesh{2,6};

%% loop over all edges aligned with x for surface contributions

for ex = 1:size(edge_topo_x,1)

    edge_nd = edge_topo_x(ex,:); % extract nodes belonging to the edge
    
    % extract only x-coordinates for the current edge
    edge_coord = coords(edge_nd,1)';
    
    % extract dof ids associated with nodes of current edge
    edge_dofs = nd_dof_map(edge_nd,:)';
    
    % compute residual and jacobian contribution from surface
    % contributions
    [res_srf, jac_srf] = compute_res_jac_srf( glb_sol(edge_dofs), edge_coord, edge_area_x(ex), ep{1} );
               
    % sum into global reidual and jacobian
    glb_res(edge_dofs)           = glb_res(edge_dofs          ) + res_srf;
    glb_jac(edge_dofs,edge_dofs) = glb_jac(edge_dofs,edge_dofs) + jac_srf;
end

%% loop over all edges aligned with y for surface contributions

for ey = 1:size(edge_topo_y,1)
    
    edge_nd = edge_topo_y(ey,:); % extract nodes belonging to the edge
    
    % extract only y-coordinates for the current edge
    edge_coord = coords(edge_nd,2)';
    
    % extract dof ids associated with nodes of current edge
    edge_dofs = nd_dof_map(edge_nd,:)';
    
    % compute residual and jacobian contribution from surface
    % contributions
    [res_srf, jac_srf] = compute_res_jac_srf( glb_sol(edge_dofs), edge_coord, edge_area_y(ey), ep{1} );
      
    % sum into global reidual and jacobian
    glb_res(edge_dofs)           = glb_res(edge_dofs          ) + res_srf;
    glb_jac(edge_dofs,edge_dofs) = glb_jac(edge_dofs,edge_dofs) + jac_srf;
end

%% loop over all nodes for volume contribution

for nd = 1:size(coords,1)
    
    % extract only coordinates for the current node
    coord = coords(nd,:);
    
    % extract dof ids associated with current node
    nd_dofs = nd_dof_map(nd,:)';
    
    % compute residual and jacobian contribution from volumetric
    % contributions
    [res_vol, jac_vol] = compute_res_jac_vol( glb_sol(nd_dofs), coord, nd_vol(nd), ep{1} );
    
    % sum into global reidual and jacobian
    glb_res(nd_dofs)         = glb_res(nd_dofs        ) + res_vol;
    glb_jac(nd_dofs,nd_dofs) = glb_jac(nd_dofs,nd_dofs) + jac_vol;
end

end