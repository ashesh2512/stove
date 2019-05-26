function [sol,fdof] = apply_bc(pp,mesh,bc,nd_dof_map)
% apply initial and boundary conditions based on problem selected
%
% Input:  ep         - element properties
%         mesh       - mesh object created uing meshgen
%         nd_dof_map - map between node and dofs
%
% Output: sol  - solution array
%         fdof - array of non-constrained degrees of freedom

% extract coordinates 
coords = mesh{2, 1};

% extract edge topologies
edge_topo_x = mesh{2,2};
edge_topo_y = mesh{2,4};

% extract sidesets
bottom_ss = mesh{2,7};
right_ss  = mesh{2,8};
top_ss    = mesh{2,9};
left_ss   = mesh{2,10};

% initialize solution array
ndof_nd = pp('dof per node');
sol = zeros(size(coords,1)*ndof_nd,1);

% collect all dirichlet nodes
dirichlet_nd = [];
if (bc('bottom') == "dirichlet")
    dirichlet_nd = [dirichlet_nd; reshape(edge_topo_x(bottom_ss,:),[],1)];
end
if (bc('right') == "dirichlet")
    dirichlet_nd = [dirichlet_nd; reshape(edge_topo_y(right_ss,:),[],1)];
end
if (bc('top') == "dirichlet")
    dirichlet_nd = [dirichlet_nd; reshape(edge_topo_x(top_ss,:),[],1)];
end
if (bc('left') == "dirichlet")
    dirichlet_nd = [dirichlet_nd; reshape(edge_topo_y(left_ss,:),[],1)];
end
dirichlet_nd = unique(dirichlet_nd); 

% apply bc based on problem
switch pp('prblm')
    case "steady heat MMS" % steady heat conduction using an MMS: T = 1/4*( cos(2*\pi*x) + cos(2*\pi*y) )
        
        cond = pp('conductivity'); % extract conductivity
                    
        sol(nd_dof_map(dirichlet_nd))= cond/4*(cos( 2*pi*coords(dirichlet_nd,1) ) ...
                                              +cos( 2*pi*coords(dirichlet_nd,2) ));

    otherwise
        error('Do not recognize the problem; check problem parameters in driver');
end

% extract array of free dofs
fdof = setdiff(nd_dof_map, nd_dof_map(dirichlet_nd));

end
