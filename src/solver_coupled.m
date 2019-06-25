function [glb_sol_np1,glb_sol_n,glb_sol_nm1] = solver_coupled(mesh_obj1,mesh_obj2,donor_map1,donor_map2,iblank1,iblank2, ...
                                                              glb_sol_np1,glb_sol_n,glb_sol_nm1,nd_dof_map1,nd_dof_map2,glb_fdof, ...
                                                              ov_info,time_info,lin_sol_info,pp)
% linear solver for coupled meshes

% extract linear solver properties
N_iters      = lin_sol_info('Newton steps');
resnrmdrop   = lin_sol_info('residual tolerance');

% perform coupled Newton solve
% loop over newton steps
for in = 1:N_iters

    % evaluate residual and jacobian matrices for all grids
    [RES, JAC_MAT] = Build_res_jac_coupled({mesh_obj1,mesh_obj2}, {donor_map1,donor_map2}, {iblank1,iblank2}, ...
                                            glb_sol_np1, glb_sol_n, glb_sol_nm1, {nd_dof_map1,nd_dof_map2}, ...
                                            ov_info, time_info, pp);

    fprintf('residual norm = %e \n', norm(RES(glb_fdof)));
    if (norm(RES(glb_fdof)) < resnrmdrop)
        break;
    elseif in == N_iters
        error('Maximum Newton iterations reached, and solution still not converged');
    end

    dsol = JAC_MAT(glb_fdof,glb_fdof)\RES(glb_fdof); % perform linear solve

    glb_sol_np1(glb_fdof) = glb_sol_np1(glb_fdof) - dsol; % update solution after linear solve

end

% update temporal solution arrays
glb_sol_nm1 = glb_sol_n; 
glb_sol_n   = glb_sol_np1;

end