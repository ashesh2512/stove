function [glb_sol_np1,glb_sol_n,glb_sol_nm1] = solver_decoupled(mesh_obj1,mesh_obj2,donor_map1,donor_map2,iblank1,iblank2, ...
                                                                glb_sol_np1,glb_sol_n,glb_sol_nm1,nd_dof_map1,nd_dof_map2,glb_fdof, ...
                                                                ov_info,time_info,lin_sol_info,pp)
% linear solver for decoupled meshes

global gmres_iter_count; % global variable to keep track of gmres iteration count
gmres_iter_count = 0; 

% extract linear solver properties
O_iters    = lin_sol_info('decoupled loops');
N_iters    = lin_sol_info('Newton steps');
if(ov_info('solve type') == "decoupled iterative" || lin_sol_info('type') == "GMRES")
    if(N_iters ~= 1)
        error('iterative solver within Newton is redundant for a decoupled scheme');
    end  
end

newton_tol = lin_sol_info('Newton tolerance');
sol_tol    = lin_sol_info('solution tolerance');

% number of points in the mesh
numpts = size(glb_sol_np1,1)/size(nd_dof_map1,2);

% loop over decoupled mesh solves
for io = 1:O_iters
    
    fprintf('decoupled loop = %d \n', io);
    
    glb_sol_np1_prev = glb_sol_np1;

    % perform Newton solve
    for in = 1:N_iters

        % evaluate residual and jacobian matrices for all grids
        [RES, JAC_MAT] = Build_res_jac({mesh_obj1,mesh_obj2}, {donor_map1,donor_map2}, {iblank1,iblank2}, ...
                                       glb_sol_np1, glb_sol_n, glb_sol_nm1, {nd_dof_map1,nd_dof_map2}, ...
                                       ov_info, time_info, pp);

        if(ov_info('solve type') ~= "decoupled iterative" && lin_sol_info('type') ~= "GMRES")
            fprintf('residual norm = %e \n', norm(RES(glb_fdof)));
            if (in == 1)
                ref_res_norm = norm(RES(glb_fdof));
                if (ref_res_norm < newton_tol)
                    break;
                end
            else
                res_norm = norm(RES(glb_fdof));
                if (res_norm/ref_res_norm < newton_tol)
                    break;
                end  
                if (in == N_iters)
                    error('Maximum Newton iterations reached, and solution still not converged');
                end
            end
        end
        
        if(ov_info('solve type') == "decoupled iterative")
            dsol = iterative_fringe_update(RES,JAC_MAT,mesh_obj1,mesh_obj2,donor_map1,donor_map2, ...
                                           glb_sol_np1,nd_dof_map1,nd_dof_map2,glb_fdof,ov_info,lin_sol_info);
        else
            if(lin_sol_info('type') == "GMRES")
                % build preconditioner for JAC_MAT
                [L,U] = ilu(JAC_MAT(glb_fdof,glb_fdof),struct('type','ilutp','droptol',1e-3));
                
                % GMRES solves
                [dsol,~,~,~,resvec] = gmres(JAC_MAT(glb_fdof,glb_fdof),RES(glb_fdof),[], ...
                                            lin_sol_info('GMRES tolerance'),lin_sol_info('GMRES iterations'),L,U); 
                                        
                fprintf('GMRES esidual norm = %e \n', resvec(end));

                iter_count = iter_count + length(resvec)-1;
                fprintf('Total GMRES iterations = %d \n\n', iter_count);
            else
               dsol = JAC_MAT(glb_fdof,glb_fdof)\RES(glb_fdof); % perform direct solve
            end
        end
        
        glb_sol_np1(glb_fdof) = glb_sol_np1(glb_fdof) - dsol; % update solution after linear solve

    end
    
    % determine solution at fringe points for all meshes
    glb_sol_np1 = fringe_interpolation({mesh_obj1,mesh_obj2}, {donor_map1,donor_map2}, ...
                                        glb_sol_np1, {nd_dof_map1,nd_dof_map2}, ov_info);
    
    L2_diff = sqrt( sum((glb_sol_np1 - glb_sol_np1_prev).^2) / numpts );
    fprintf('solution diff = %e \n\n', L2_diff);
    
    if (L2_diff < sol_tol) % check if glb_sol_np1 is still changing
        break;
    end

end

% update temporal solution arrays
glb_sol_nm1 = glb_sol_n; 
glb_sol_n   = glb_sol_np1;

end


%% decoupled GMRES solve with intermediate fringe updates

function dsol = iterative_fringe_update(RES,JAC_MAT,mesh_obj1,mesh_obj2,donor_map1,donor_map2, ...
                                        fulldsol,nd_dof_map1,nd_dof_map2,glb_fdof,ov_info,lin_sol_info)
global gmres_iter_count; 

gmres_tol   = lin_sol_info('GMRES tolerance');
maxit       = lin_sol_info('GMRES iterations');
dsol_tol    = lin_sol_info('dsol tolerance');
max_exch_it = ov_info('exchange iterations');

numpts = size(glb_fdof,1);
dsol = zeros(numpts,1); % initialize dsol array

% build preconditioner for JAC_MAT
[L,U] = ilu(JAC_MAT(glb_fdof,glb_fdof),struct('type','ilutp','droptol',1e-3));

loc_count = 0; % initialize iteration counter

fprintf('\n');

while (loc_count <= max_exch_it)
    
    fulldsol(glb_fdof,1) = dsol; % update full delta solution array
    dsol_prev = dsol; % update previous dsol array
    
    % determine RHS correction at fringe points for all meshes
    RES = fringe_DBC(RES, {mesh_obj1,mesh_obj2}, {donor_map1,donor_map2}, ...
                     fulldsol, {nd_dof_map1,nd_dof_map2}, ov_info);
    
    % GMRES solves
    [dsol,~,~,~,resvec] = gmres(JAC_MAT(glb_fdof,glb_fdof),RES(glb_fdof),[], ...
                                gmres_tol,maxit,L,U,fulldsol(glb_fdof));                                                      
    loc_count = loc_count+(length(resvec)-1);
    if(resvec(end) < eps) % break if absolute residual is very low
        break;
    end

    dsol_diff = sqrt( sum((dsol - dsol_prev).^2) / numpts );
    if(dsol_diff < dsol_tol) % check if dsol is still changing
        break
    end

    fprintf('GMRES esidual norm = %e \n', resvec(end));
    
end

gmres_iter_count = gmres_iter_count + loc_count;
fprintf('Total GMRES iterations = %d \n\n', gmres_iter_count);

end
 



