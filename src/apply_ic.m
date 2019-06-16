function sol = apply_ic(pp,mesh,nd_dof_map,init_time)
% apply initial and boundary conditions based on problem selected
%
% Input:  pp         - element properties
%         mesh       - mesh object created uing meshgen
%         nd_dof_map - map between node and dofs
%         init_time  - initial time

% Output: sol  - solution array

% extract coordinates 
coords = mesh{2, 1};

% initialize solution array
ndof_nd = pp('dof per node');
sol = zeros(size(coords,1)*ndof_nd,1);

% apply bc based on problem
switch pp('prblm')
    case "steady heat MMS"
        return;
        
    case "unsteady scalar adv" % sin(freq*X - vel*time) + sin(freq*Y - vel*time);
        
        freq = pp('frequency'); % extract frequency
        vel  = pp('velocity');  % extract constant velocity
        
        sol = sin(freq*coords(:,1) - vel*init_time) + sin(freq*coords(:,2) - vel*init_time);
        
    otherwise
        error('Do not recognize the problem; check problem parameters in driver');
end

end
