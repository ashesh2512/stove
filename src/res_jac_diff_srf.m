function [res,jac] = res_jac_diff_srf(sol, edge_coord, edge_area, pp)
% compute surface residual and jacobian contribution for steady heat
% diffusion

% initialize arrays
dofs = length(sol);
res  = zeros(dofs,1);
jac  = zeros(dofs,dofs);

% conductivity
cond = pp('conductivity');

% q = -kdt/dx
q      = -cond*(sol(2) - sol(1))/abs(edge_coord(2) - edge_coord(1));
qcoeff = -cond                  /abs(edge_coord(2) - edge_coord(1));
 
% assign resiual contribution => q.ndS
res(1) = res(1) + q*(+edge_area); % left node
res(2) = res(2) + q*(-edge_area); % right node

% assemble jacobian
jac(1,1) = jac(1,1) + -qcoeff*(+edge_area); % left  w.r.t. left
jac(1,2) = jac(1,2) +  qcoeff*(+edge_area); % left  w.r.t. right
jac(2,1) = jac(2,1) + -qcoeff*(-edge_area); % right w.r.t. left
jac(2,2) = jac(2,2) +  qcoeff*(-edge_area); % right w.r.t. right

end