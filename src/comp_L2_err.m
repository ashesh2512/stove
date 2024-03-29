function L2_err = comp_L2_err(pp,sol,mesh1,mesh2,time)
% This function computes the L2 error and plots the analytical and
% discretized solutions for qualitative comparison

% initialize L2 error array
num_dof_nd = pp('dof per node');
L2_err = zeros(num_dof_nd,1);

% number of points in the mesh
coords1 = mesh1{2,1};
coords2 = mesh2{2,1};
coords  = [coords1;coords2];
numpts  = size(coords,1);

switch pp('prblm')
    case "Laplace"
        
        sol_anlyt = sin(coords(:,1)).*sinh(coords(:,2));
        
        L2_err(1) = sqrt( sum((sol - sol_anlyt).^2) / numpts );
        
    case "steady heat MMS" % steady heat conduction using an MMS: T = k/4(cos(2*\pi*x) + cos(2*\pi*y))
        
        cond = pp('conductivity');
        
        sol_anlyt = cond/4*(cos(2*pi*coords(:,1)) + cos(2*pi*coords(:,2)));
        
        L2_err(1) = sqrt( sum((sol - sol_anlyt).^2) / numpts );
        
    case "unsteady scalar adv" % sin(freq*X - vel*time) + sin(freq*Y - vel*time);
        % analytical solution is different than for "unsteady adv diff"
        % because the inertial term is different - FIXME

        vel = pp('velocity');

        sol_anlyt = sin(pi*coords(:,1) - vel(1)*time) ...
                  + sin(pi*coords(:,2) - vel(2)*time);
        
        L2_err(1) = sqrt( sum((sol - sol_anlyt).^2) / numpts );
        
    case "unsteady adv diff"
        
        cond = pp('conductivity'); % extract conductivity
        vel  = pp('velocity');  % extract constant velocity
        avg  = pp('mean'); % extr
        
        sol_anlyt = exp(-cond*(pi^2)*time)...
            *(sin(pi*(coords(:,1)-vel(1)*time)) + sin(pi*(coords(:,2)-vel(2)*time)))...
            + avg;
        
        L2_err(1) = sqrt( sum((sol - sol_anlyt).^2) / numpts );
        
    otherwise
        error('Do not recognize the problem id');
end

fprintf('\nL2 error = %16.16e', L2_err);
fprintf('\n');

end
