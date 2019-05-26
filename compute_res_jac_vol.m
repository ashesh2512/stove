function [res_vol,jac_vol] = compute_res_jac_vol(sol, coord, nd_vol, pp)
% This function computes the volumetric conribution to the residual and jacobain.

% initialize arrays
numdofs = size(sol,1); % total number of dofs
res_vol = zeros(numdofs,1);       % initialize residual controbution
jac_vol = zeros(numdofs,numdofs); % initialize jacobian controbution

switch pp('prblm')
    case "steady heat MMS" % steady heat conduction using an MMS: q^dot = k*pi^2*(cos(2*\pi*x) + cos(2*\pi*y))
        
        cond = pp('conductivity'); % conductivity
        
        x = coord(1);
        y = coord(2);
        
        res_vol = res_vol + cond*(pi^2)*(cos(2.0*pi*x) + cos(2.0*pi*y))*nd_vol;
        
    otherwise
        error('Do not recognize the problem id');
end

end
