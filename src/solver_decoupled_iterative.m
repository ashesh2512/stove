function [glb_sol_np1,glb_sol_n,glb_sol_nm1] = solver_decoupled_iterative(mesh_obj1,mesh_obj2,donor_map1,donor_map2,iblank1,iblank2, ...
                                                                          glb_sol_np1,glb_sol_n,glb_sol_nm1,nd_dof_map1,nd_dof_map2,glb_fdof,glb_fld_dof, ...
                                                                          ov_info,time_info,lin_sol_info,pp)
% linear solver for decoupled meshes

global gmres_iter_count; % global variable to keep track of gmres iteration count
gmres_iter_count = 0; 

% extract linear solver properties
O_iters = lin_sol_info('decoupled loops');
sol_tol = lin_sol_info('solution tolerance');

% number of points in the mesh
numpts = size(glb_sol_np1,1)/size(nd_dof_map1,2);

% loop over decoupled mesh solves
for io = 1:O_iters
    
    fprintf('decoupled loop = %d \n', io);
    
    glb_sol_np1_prev = glb_sol_np1; % update previous solution array

    % evaluate residual and jacobian matrices for all grids
    [RES, JAC_MAT] = Build_res_jac({mesh_obj1,mesh_obj2}, {donor_map1,donor_map2}, {iblank1,iblank2}, ...
                                   glb_sol_np1, glb_sol_n, glb_sol_nm1, {nd_dof_map1,nd_dof_map2}, ...
                                   ov_info, time_info, pp);

    if(ov_info('fringe update') == "iterative") %
        fulldsol = iterative_fringe_update(RES,JAC_MAT,mesh_obj1,mesh_obj2,donor_map1,donor_map2, ...
                                           nd_dof_map1,nd_dof_map2,glb_fdof,ov_info,lin_sol_info);
        
        glb_sol_np1(glb_fld_dof) = glb_sol_np1(glb_fld_dof) - fulldsol(glb_fld_dof); % update solution after linear solve
    else
        % build preconditioner for JAC_MAT
        [L,U] = ilu(JAC_MAT(glb_fdof,glb_fdof),struct('type','ilutp','droptol',1e-3));
        
        % GMRES solves
        [dsol,~,~,~,resvec] = gmres(JAC_MAT(glb_fdof,glb_fdof),RES(glb_fdof),[], ...
                                    lin_sol_info('GMRES tolerance'),lin_sol_info('GMRES iterations'),L,U);
        
        fprintf('GMRES esidual norm = %e \n', resvec(end));
        
        gmres_iter_count = gmres_iter_count + length(resvec)-1;
        fprintf('Total GMRES iterations = %d \n\n', gmres_iter_count);
        
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

function fulldsol = iterative_fringe_update(RES,JAC_MAT,mesh_obj1,mesh_obj2,donor_map1,donor_map2, ...
                                            nd_dof_map1,nd_dof_map2,glb_fdof,ov_info,lin_sol_info)
global gmres_iter_count; 

gmres_tol   = lin_sol_info('GMRES tolerance');
maxit       = lin_sol_info('GMRES iterations');
dsol_tol    = lin_sol_info('dsol tolerance');
max_exch_it = ov_info('exchange iterations');

% initialize dull change in solution array 
fulldsol = zeros(size(RES,1));
dsol = fulldsol(glb_fdof);

% build preconditioner for JAC_MAT
[L,U] = ilu(JAC_MAT(glb_fdof,glb_fdof),struct('type','ilutp','droptol',1e-3));

loc_count = 0; % initialize iteration counter

fprintf('\n');

while (loc_count < max_exch_it)

    dsol_prev = dsol; % update previous dsol array
    
    % determine RHS correction at fringe points for all meshes
    RES = fringe_DBC(RES, {mesh_obj1,mesh_obj2}, {donor_map1,donor_map2}, ...
                     fulldsol, {nd_dof_map1,nd_dof_map2}, ov_info);
    
    % GMRES solves
    [dsol,~,~,~,resvec] = gmres(JAC_MAT(glb_fdof,glb_fdof),RES(glb_fdof),[], ...
                                gmres_tol,maxit,L,U,fulldsol(glb_fdof,1));
    loc_count = loc_count+(length(resvec)-1);
    
    fulldsol(glb_fdof) = dsol; % update change in solution array
    
    fprintf('GMRES esidual norm = %e \n', resvec(end));    
    if(resvec(end) < eps) % break if absolute residual is very low
        break;
    end

    dsol_diff = sqrt( sum((dsol - dsol_prev).^2) / size(glb_fdof,1) );
    fprintf('dsol diff = %e \n', dsol_diff);
    if(dsol_diff < dsol_tol) % check if dsol is still changing
        break
    end
    
end

gmres_iter_count = gmres_iter_count + loc_count;
fprintf('Total GMRES iterations = %d \n\n', gmres_iter_count);

end
 



