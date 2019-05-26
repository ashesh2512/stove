function [res, jac] = compute_res_jac_srf(sol, edge_coord, edge_area, ep)
% This function evaluates the scs contribution to the nodal residual and jacobian.

% determine case
prblm = ep(1);

switch prblm
    case 1 % steady heat conduction
        
        [res,jac] = res_jac_diff_srf(sol, edge_coord, edge_area, ep);
        
    otherwise
        error('Do not recognize the problem id');
end 

end