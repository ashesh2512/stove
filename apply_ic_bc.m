function [sol,fdof] = apply_ic_bc(ep,mesh,nd_dof_map)
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
ndof_nd = ep(2);
sol = zeros(size(coords,1)*ndof_nd,1);

% problem id
prblm = ep(1);

switch prblm
    case 1 % steady heat conduction using an MMS: T = 1/4*( cos(2*\pi*x) + cos(2*\pi*y) )
        
        cond = ep(3); % extract conductivity
        
        % dirichlet bc applied to all four boundaries
        dirichlet_nd = [reshape(edge_topo_x(bottom_ss,:),[],1); ...
                        reshape(edge_topo_y(right_ss,:),[],1); ...
                        reshape(edge_topo_x(top_ss,:),[],1); ...
                        reshape(edge_topo_y(left_ss,:),[],1)];
        dirichlet_nd = unique(dirichlet_nd); 
                    
        sol(nd_dof_map(dirichlet_nd))= cond/4*(cos( 2*pi*coords(dirichlet_nd,1) ) ...
                                              +cos( 2*pi*coords(dirichlet_nd,2) ));

        
        fdof = setdiff(nd_dof_map, nd_dof_map(dirichlet_nd));
        
    otherwise
        error('Do not recognize the problem id');
end

end
