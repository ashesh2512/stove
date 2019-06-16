function sol = apply_bc(pp,mesh,cnd,nd_dof_map,time)
% apply initial and boundary conditions based on problem selected
%
% Input:  pp         - element properties
%         mesh       - mesh object created uing meshgen
%         cnd        - lit of constrained nodes
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
    case "steady heat MMS" % steady heat conduction using an MMS: T = 1/4*( cos(2*\pi*x) + cos(2*\pi*y) )
        
        cond = pp('conductivity'); % extract conductivity
                    
        sol(nd_dof_map(cnd))= cond/4*(cos( 2*pi*coords(cnd,1) ) ...
                                     +cos( 2*pi*coords(cnd,2) ));

    case "unsteady scalar adv" % sin(freq*X - vel*time) + sin(freq*Y - vel*time);
        
        freq = pp('frequency'); % extract frequency
        vel  = pp('velocity');  % extract constant velocity
        
        sol(nd_dof_map(cnd)) = sin(freq*coords(cnd,1) - vel*time) + sin(freq*coords(cnd,2) - vel*time);
        
    otherwise
        error('Do not recognize the problem; check problem parameters in driver');
end

end