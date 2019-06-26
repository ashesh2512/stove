function [glb_sol_np1,glb_sol_n,glb_sol_nm1] = solver_decoupled(mesh_obj1,mesh_obj2,donor_map1,donor_map2,iblank1,iblank2, ...
                                                                glb_sol_np1,glb_sol_n,glb_sol_nm1,nd_dof_map1,nd_dof_map2,glb_fdof, ...
                                                                ov_info,time_info,lin_sol_info,pp)
% linear solver for decoupled meshes

% extract linear solver properties
O_iters    = lin_sol_info('decoupled loops');
N_iters    = lin_sol_info('Newton steps');
resnrmdrop = lin_sol_info('residual tolerance');

% loop over decoupled mesh solves
for io = 1:O_iters
    
    fprintf('decoupled loop = %d \n', io);

    % perform Newton solve
    for in = 1:N_iters

        % evaluate residual and jacobian matrices for all grids
        [RES, JAC_MAT] = Build_res_jac({mesh_obj1,mesh_obj2}, {donor_map1,donor_map2}, {iblank1,iblank2}, ...
                                       glb_sol_np1, glb_sol_n, glb_sol_nm1, {nd_dof_map1,nd_dof_map2}, ...
                                       ov_info, time_info, pp);

        fprintf('residual norm = %e \n', norm(RES(glb_fdof)));
        if (in == 1)
            ref_res_norm = norm(RES(glb_fdof));
            if (ref_res_norm < resnrmdrop)
                break;
            end
        else
            res_norm = norm(RES(glb_fdof));
            if (res_norm/ref_res_norm < resnrmdrop)
                break;
            end  
            if (in == N_iters)
                error('Maximum Newton iterations reached, and solution still not converged');
            end
        end

        dsol = JAC_MAT(glb_fdof,glb_fdof)\RES(glb_fdof); % perform linear solve

        glb_sol_np1(glb_fdof) = glb_sol_np1(glb_fdof) - dsol; % update solution after linear solve

    end
    
    % determine solution at fringe points for all meshes
    glb_sol_np1 = fringe_interpolation({mesh_obj1,mesh_obj2}, {donor_map1,donor_map2}, ...
                                        glb_sol_np1, {nd_dof_map1,nd_dof_map2}, ov_info);
    
    fprintf('\n');

end

% update temporal solution arrays
glb_sol_nm1 = glb_sol_n; 
glb_sol_n   = glb_sol_np1;

end


