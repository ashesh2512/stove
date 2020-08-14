function sol = apply_ic_n(pp,mesh,init_time)
% apply initial and boundary conditions based on problem selected
%
% Input:  pp         - element properties
%         mesh       - mesh object created uing meshgen
%         init_time  - initial time

% Output: sol  - solution array

% extract coordinates 
coords = mesh{2, 1};

% initialize solution array
ndof_nd = pp('dof per node');
sol = zeros(size(coords,1)*ndof_nd,1);

% apply bc based on problem
switch pp('prblm')
    case {"steady heat MMS","Laplace"}
        return;
        
    case "unsteady scalar adv" % sin(pi*X - vel*time) + sin(pi*Y - vel*time)
        
        vel  = pp('velocity');  % extract constant velocity
        
        sol = sin(pi*coords(:,1) - vel(1)*init_time) + sin(pi*coords(:,2) - vel(2)*init_time);
        
    case "unsteady adv diff"
        
        cond = pp('conductivity'); % extract conductivity
        vel  = pp('velocity');  % extract constant velocity
        avg  = pp('mean'); % extr
        
        sol = exp(-cond*(pi^2)*init_time)...
            *(sin(pi*(coords(:,1)-vel(1)*init_time)) + sin(pi*(coords(:,2)-vel(2)*init_time)))...
            + avg;
        
    otherwise
        error('Do not recognize the problem; check problem parameters in driver');
end

end
