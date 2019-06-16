function [res,jac] = res_jac_scalar_adv_srf(sol,edge_coord,edge_area,edge_nrml,pp)
% compute surface residual and jacobian contribution for scalar advection

% initialize arrays
dofs = length(sol);
res  = zeros(dofs,1);
jac  = zeros(dofs,dofs);

% extract advection constants
vel = pp('velocity');
vel = dot(vel,edge_nrml);

% interpolate transport variable at edge midpoint
avg_sol = (sol(1) + sol(2))/2;
 
% assign resiual contribution => \phi*v.ndS
res(1) = res(1) + vel*avg_sol*(+edge_area); % left node
res(2) = res(2) + vel*avg_sol*(-edge_area); % right node

% assemble jacobian
jac(1,1) = jac(1,1) + vel*0.5*(+edge_area); % left  w.r.t. left
jac(1,2) = jac(1,2) + vel*0.5*(+edge_area); % left  w.r.t. right
jac(2,1) = jac(2,1) + vel*0.5*(-edge_area); % right w.r.t. left
jac(2,2) = jac(2,2) + vel*0.5*(-edge_area); % right w.r.t. right

end