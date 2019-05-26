function [res_vol,jac_vol] = compute_res_jac_vol(sol, coord, nd_vol, ep)
% This function computes the volumetric conribution to the residual and jacobain.

% initialize arrays
numdofs = size(sol,1); % total number of dofs
res_vol = zeros(numdofs,1);       % initialize residual controbution
jac_vol = zeros(numdofs,numdofs); % initialize jacobian controbution

% determine case
prblm   = ep(1); % problem type

switch prblm
    case 1 % steady heat conduction using an MMS: q^dot = k*pi^2*(cos(2*\pi*x) + cos(2*\pi*y))
        
        cond = ep(3); % conductivity
        
        x = coord(1);
        y = coord(2);
        
        res_vol = res_vol + cond*(pi^2)*(cos(2.0*pi*x) + cos(2.0*pi*y))*nd_vol;
        
    otherwise
        error('Do not recognize the problem id');
end

end
