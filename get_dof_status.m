function [fdof,const_nd] = get_dof_status(mesh,bc,nd_dof_map)
% apply initial and boundary conditions based on problem selected
%
% Input:  mesh       - mesh object created uing meshgen
%         bc         - boundary conditions
%         nd_dof_map - map between node and dofs
%
% Output: sol      - solution array
%         fdof     - array of non-constrained degrees of freedom
%         const_nd - array of constrained nodes

% extract edge topologies
edge_topo_x = mesh{2,2};
edge_topo_y = mesh{2,4};

% extract sidesets
bottom_ss = mesh{2,7};
right_ss  = mesh{2,8};
top_ss    = mesh{2,9};
left_ss   = mesh{2,10};

% collect all constrained nodes
const_nd = [];
if (bc('bottom') == "dirichlet")
    const_nd = [const_nd; reshape(edge_topo_x(bottom_ss,:),[],1)];
end
if (bc('right') == "dirichlet")
    const_nd = [const_nd; reshape(edge_topo_y(right_ss,:),[],1)];
end
if (bc('top') == "dirichlet")
    const_nd = [const_nd; reshape(edge_topo_x(top_ss,:),[],1)];
end
if (bc('left') == "dirichlet")
    const_nd = [const_nd; reshape(edge_topo_y(left_ss,:),[],1)];
end
const_nd = unique(const_nd); 

% extract array of free dofs
fdof = setdiff(nd_dof_map, nd_dof_map(const_nd));

end
