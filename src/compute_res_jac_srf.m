function [res, jac] = compute_res_jac_srf(sol, edge_coord, edge_area, edge_nrml, pp)
% This function evaluates the scs contribution to the nodal residual and jacobian.

switch pp('prblm')
    case "steady heat MMS" % steady heat conduction
        
        [res,jac] = res_jac_diff_srf(sol, edge_coord, edge_area, pp);
        
    case "unsteady scalar adv" % unsteady scalar advection with constant velocity
        
        [res,jac] = res_jac_scalar_adv_srf(sol, edge_coord, edge_area, edge_nrml, pp);

    otherwise
        error('Do not recognize the problem id');
end 

end