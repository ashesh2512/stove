function sol = apply_dbc(pp,mesh,dnd,nd_dof_map,time)
% apply initial and boundary conditions based on problem selected
%
% Input:  pp         - element properties
%         mesh       - mesh object created uing meshgen
%         dnd        - lit of constrained dirichlet nodes
%         nd_dof_map - map between node and dofs
%
% Output: sol  - solution array

% extract coordinates 
coords = mesh{2, 1};

% initialize solution array
ndof_nd = pp('dof per node');
sol = zeros(size(coords,1)*ndof_nd,1);

% apply bc based on problem
switch pp('prblm')
    case "Laplace"
        sol(nd_dof_map(dnd)) = sin(coords(dnd,1)).*sinh(coords(dnd,2));
        
    case "steady heat MMS" % steady heat conduction using an MMS: T = 1/4*( cos(2*\pi*x) + cos(2*\pi*y) )
        
        cond = pp('conductivity'); % extract conductivity
                    
        sol(nd_dof_map(dnd))= cond/4*(cos( 2*pi*coords(dnd,1) ) ...
                                     +cos( 2*pi*coords(dnd,2) ));

    case "unsteady scalar adv" % sin(pi*X - vel*time) + sin(pi*Y - vel*time)
                
        vel  = pp('velocity');  % extract constant velocity
        
        % If dnd is empty, which is the case for pnd, solution array
        % continues to be consistent
        sol(nd_dof_map(dnd)) = sin(pi*coords(dnd,1) - vel(1)*time) ...
                             + sin(pi*coords(dnd,2) - vel(2)*time);
                         
    case "unsteady adv diff"
        if(~isempty(dnd))
            error('Dirichlet BC not setup for unsteady advection diffusion case');
        end
        
    otherwise
        error('Do not recognize the problem; check problem parameters in driver');
end

end
