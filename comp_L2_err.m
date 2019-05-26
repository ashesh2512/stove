function L2_err = comp_L2_err(ep,sol,coord)
% This function computes the L2 error and plots the analytical and
% discretized solutions for qualitative comparison

% initialize L2 error array
num_dof_nd = ep(2);
L2_err = zeros(num_dof_nd,1);

% number of points in the mesh
numpts = size(coord,1);

% problem type
prblm = ep(1);

switch prblm
    case 1 % steady heat conduction using an MMS: T = k/4(cos(2*\pi*x) + cos(2*\pi*y))
        
        cond = ep(3);
        
        xsol_anlyt = cond/4*(cos(2*pi*coord(:,1)) + cos(2*pi*coord(:,2)));
        
        L2_err(1) = sqrt( sum((sol - xsol_anlyt).^2) / numpts );

    otherwise
        error('Do not recognize the problem id');
end

end
