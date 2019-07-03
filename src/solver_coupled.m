function [glb_sol_np1,glb_sol_n,glb_sol_nm1] = solver_coupled(mesh_obj1,mesh_obj2,donor_map1,donor_map2,iblank1,iblank2, ...
                                                              glb_sol_np1,glb_sol_n,glb_sol_nm1,nd_dof_map1,nd_dof_map2,glb_fdof, ...
                                                              ov_info,time_info,lin_sol_info,pp)
% linear solver for coupled meshes

% extract linear solver properties
N_iters      = lin_sol_info('Newton steps');
resnrmdrop   = lin_sol_info('Newton tolerance');

gmres_iter_count = 0; % initialize gmres iteration counter

% perform coupled Newton solve
% loop over newton steps
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

    if(lin_sol_info('type') == "GMRES")
        % build preconditioner for JAC_MAT
        [L,U] = ilu(JAC_MAT(glb_fdof,glb_fdof),struct('type','ilutp','droptol',1e-3));

        % GMRES solves
        [dsol,~,~,~,resvec] = gmres(JAC_MAT(glb_fdof,glb_fdof),RES(glb_fdof),[], ...
                                    lin_sol_info('GMRES tolerance'),lin_sol_info('GMRES iterations'),L,U); 
        gmres_iter_count = gmres_iter_count + length(resvec)-1;
        
        fprintf('Total GMRES iterations = %d \n\n', gmres_iter_count);
    else
       dsol = JAC_MAT(glb_fdof,glb_fdof)\RES(glb_fdof); % perform direct solve
    end
            
    glb_sol_np1(glb_fdof) = glb_sol_np1(glb_fdof) - dsol; % update solution after linear solve

end

% update temporal solution arrays
glb_sol_nm1 = glb_sol_n; 
glb_sol_n   = glb_sol_np1;

end