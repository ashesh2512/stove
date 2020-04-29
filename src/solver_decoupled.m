function [glb_sol_np1,glb_sol_n,glb_sol_nm1] = solver_decoupled(mesh_obj1,mesh_obj2,mesh_obj1_cell,mesh_obj2_cell, ...
                                                                donor_map1,donor_map2,iblank1,iblank2,iblank_cell1,iblank_cell2, ...
                                                                glb_sol_np1,glb_sol_n,glb_sol_nm1,nd_dof_map1,nd_dof_map2,cell_dof_map1,cell_dof_map2,glb_fdof, ...
                                                                ov_info,time_info,lin_sol_info,pp)
% linear solver for decoupled meshes
L2_err = [];

% extract linear solver properties
O_iters    = lin_sol_info('decoupled loops');
N_iters    = lin_sol_info('Newton steps');
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
%             if (in == N_iters)
%                 error('Maximum Newton iterations reached, and solution still not converged');
%             end
        end

        dsol = JAC_MAT(glb_fdof,glb_fdof)\RES(glb_fdof); % perform direct solve

        glb_sol_np1(glb_fdof) = glb_sol_np1(glb_fdof) - dsol; % update solution after linear solve
        
    end
    
    % determine solution at fringe points for all meshes
    if (ov_info('fringe node averaging'))
%         coords1 = mesh_obj1{2,1};
%         glb_sol_np1(nd_dof_map1,1) = coords1(:,1) + coords1(:,2);
%         coords2 = mesh_obj2{2,1};
%         glb_sol_np1(nd_dof_map2,1) = coords2(:,1) + coords2(:,2);
%         plot_sol(mesh_obj1,mesh_obj2,glb_sol_np1,nd_dof_map1,nd_dof_map2,0.0);
        
        glb_sol_np1_cell = zeros(max([cell_dof_map1;cell_dof_map2],[],'all'),1);
        glb_sol_np1_cell = shift_solution({mesh_obj1,mesh_obj2},{iblank1,iblank2},glb_sol_np1,glb_sol_np1_cell, ...
                                          {nd_dof_map1,nd_dof_map2},{cell_dof_map1,cell_dof_map2},"cell",ov_info);
%         plot_sol_cell(mesh_obj1,glb_sol_np1_cell,cell_dof_map1,0.0);
        
        % NOTE: size of solution array corresponds to both meshes
        if(ov_info('background grid') == 1) 
            glb_sol_np1_cell = fringe_interpolation_mixed(mesh_obj1_cell,mesh_obj2,donor_map1,glb_sol_np1_cell,glb_sol_np1,cell_dof_map1,nd_dof_map2,ov_info);
            glb_sol_np1 = fringe_interpolation_mixed(mesh_obj2,mesh_obj1_cell,donor_map2,glb_sol_np1,glb_sol_np1_cell,nd_dof_map2,cell_dof_map1,ov_info);
        else
            glb_sol_np1_cell = fringe_interpolation_mixed(mesh_obj2_cell,mesh_obj1,donor_map2,glb_sol_np1_cell,glb_sol_np1,cell_dof_map2,nd_dof_map1,ov_info);
            glb_sol_np1 = fringe_interpolation_mixed(mesh_obj1,mesh_obj2_cell,donor_map1,glb_sol_np1,glb_sol_np1_cell,nd_dof_map1,cell_dof_map2,ov_info);
        end

        glb_sol_np1 = shift_solution({mesh_obj1,mesh_obj2},{iblank1,iblank2},glb_sol_np1,glb_sol_np1_cell, ...
                                     {nd_dof_map1,nd_dof_map2},{cell_dof_map1,cell_dof_map2},"node",ov_info);
%         plot_sol(mesh_obj1,mesh_obj2,glb_sol_np1,nd_dof_map1,nd_dof_map2,0.0);

%         diff_norm = norm(glb_sol_np1_avgd(nd_dof_map1,1) - glb_sol_np1(nd_dof_map1,1))
    else
        glb_sol_np1 = fringe_interpolation({mesh_obj1,mesh_obj2}, {donor_map1,donor_map2}, ...
                                           glb_sol_np1, {nd_dof_map1,nd_dof_map2}, ov_info);
    end
    
    L2_diff = sqrt( sum((glb_sol_np1 - glb_sol_np1_prev).^2) / numpts );
    fprintf('solution diff = %e \n\n', L2_diff);
    
    L2_err = [L2_err, comp_L2_err(pp,glb_sol_np1,mesh_obj1,mesh_obj2,0.0)];

    if (L2_diff < sol_tol) % check if glb_sol_np1 is still changing
        break;
    end

end

% update temporal solution arrays
glb_sol_nm1 = glb_sol_n; 
glb_sol_n   = glb_sol_np1;

end