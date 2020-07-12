function [res_vol,jac_vol] = compute_res_jac_vol(sol_np1,sol_n,sol_nm1,coord,nd_vol,time_info,pp)
% This function computes the volumetric conribution to the residual and jacobain.

% initialize arrays
numdofs = size(sol_np1,1);        % total number of dofs
res_vol = zeros(numdofs,1);       % initialize residual controbution
jac_vol = zeros(numdofs,numdofs); % initialize jacobian controbution

% extract time info
g_np1 = time_info('gamma_np1');
g_n   = time_info('gamma_n1');
g_nm1 = time_info('gamma_nm1');
dt    = time_info('time step');

switch pp('prblm')
    case "steady heat MMS" % steady heat conduction using an MMS: q^dot = k*pi^2*(cos(2*\pi*x) + cos(2*\pi*y))
        
        cond = pp('conductivity'); % conductivity
        
        x = coord(1);
        y = coord(2);
        
        res_vol(1) = res_vol - cond*(pi^2)*(cos(2.0*pi*x) + cos(2.0*pi*y))*nd_vol;
        
    case "unsteady scalar adv" % unsteady scalar advection with constant velocity
                
        res_vol(1) = res_vol(1) + pi*(g_np1*sol_np1(1) + g_n*sol_n(1) + g_nm1*sol_nm1(1))*nd_vol/dt;

        jac_vol(1,1) = jac_vol(1,1) + pi*g_np1*nd_vol/dt;
        
    case "unsteady adv diff" % unsteady scalar advection with constant velocity
                
        res_vol(1) = res_vol(1) + (g_np1*sol_np1(1) + g_n*sol_n(1) + g_nm1*sol_nm1(1))*nd_vol/dt;

        jac_vol(1,1) = jac_vol(1,1) + g_np1*nd_vol/dt;

    otherwise
        error('Do not recognize the problem id');
end

end
