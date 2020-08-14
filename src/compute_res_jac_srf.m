function [res, jac] = compute_res_jac_srf(sol, edge_coord, edge_area, edge_nrml, pp)
% This function evaluates the scs contribution to the nodal residual and jacobian.

switch pp('prblm')
    case {"steady heat MMS","Laplace"} % steady heat conduction
        
        [res,jac] = res_jac_diff_srf(sol, edge_coord, edge_area, pp);
        
    case "unsteady scalar adv" % unsteady scalar advection with constant velocity
        
        [res,jac] = res_jac_scalar_adv_srf(sol, edge_coord, edge_area, edge_nrml, pp);
        
    case "unsteady adv diff" % unsteady scalar advection diffusion with constant velocity
        
        [res_diff,jac_diff] = res_jac_diff_srf(sol, edge_coord, edge_area, pp);
        [res_adv , jac_adv] = res_jac_scalar_adv_srf(sol, edge_coord, edge_area, edge_nrml, pp);
        
        res = res_diff + res_adv;
        jac = jac_diff + jac_adv;

    otherwise
        error('Do not recognize the problem id');
end 

end