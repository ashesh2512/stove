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
                                           glb_sol_np1,nd_dof_map1,nd_dof_map2,glb_fdof,ov_info,lin_sol_info,glb_fld_dof,pp);
        
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

function fulldsol = iterative_fringe_update(RHS,JAC_MAT,mesh_obj1,mesh_obj2,donor_map1,donor_map2, ...
                                            glb_sol_np1,nd_dof_map1,nd_dof_map2,glb_fdof,ov_info,lin_sol_info,glb_fld_dof,pp)
global gmres_iter_count; 

gmres_tol   = lin_sol_info('GMRES tolerance');
maxit       = lin_sol_info('GMRES iterations');
U_dim       = lin_sol_info('dsol history');
dsol_tol    = lin_sol_info('dsol tolerance');
max_exch_it = ov_info('exchange iterations');

% initialize full change in solution array 
fulldsol = zeros(size(glb_sol_np1,1));
free_dsol = fulldsol(glb_fdof);

Q = zeros(size(glb_fdof,1), U_dim); % Orthonormal basis for the RHS subspace
U = zeros(size(glb_fdof,1), U_dim); % Solutions satisfying A*U = Q
R = zeros(           U_dim, U_dim); % R factor in the QR decomposition for the RHS subspace
T = zeros(           U_dim, U_dim); % R factor in the QR decomposition for the RHS subspace

% build preconditioner for JAC_MAT
free_JAC_MAT = JAC_MAT(glb_fdof,glb_fdof);
[L_pre,U_pre] = ilu(free_JAC_MAT,struct('type','ilutp','droptol',1e-3));

U_d = 0; % initialize counter for solution history
loc_count = 0; % initialize GMRES iteration counter
exch_it   = 0; % initialize iteration counter

glb_sol_tmp = glb_sol_np1;

L2 = [];

fprintf('\n');

while (exch_it < max_exch_it)

    dsol_prev = free_dsol; % update previous dsol array
    
    % determine RHS correction at fringe points for all meshes
    RHS_adj = fringe_DBC(RHS, JAC_MAT, {mesh_obj1,mesh_obj2}, {donor_map1,donor_map2}, ...
                         fulldsol, {nd_dof_map1,nd_dof_map2}, ov_info);
    free_RES = RHS_adj(glb_fdof);
    
    if(U_d == 0)    
        % GMRES solve
        [free_dsol,~,~,~,resvec] = gmres(free_JAC_MAT,free_RES,[],gmres_tol,maxit,L_pre,U_pre,free_dsol);
        
        orthg_RES = free_JAC_MAT*free_dsol;
        
        R(1, 1) = norm(orthg_RES);
        Q(:, 1) = orthg_RES/R(1, 1);
        U(:, 1) = free_dsol/R(1, 1);
        
    else
        % Project onto orthogonal complement of RHSs subspace.
        alpha = R(1:U_d,U_d);
        RES_tilde = free_RES - Q(:,1:U_d)*alpha;
        
        % Component of initial solution lying in RHS subspace.
        free_dsol_par = U(:,1:U_d)*alpha;
        
        % Components of initial solution perpendicular to RHS subspace.
        free_dsol_per = U_pre\(L_pre\RES_tilde);
        
        init_dsol = free_dsol_par + free_dsol_per;
        init_dsol = init_dsol + U_pre\(L_pre\(free_RES - free_JAC_MAT*init_dsol));
        
        % apply Hegedus scaling
        xi_min = free_RES'*(free_JAC_MAT*init_dsol)/(norm(free_JAC_MAT*init_dsol)^2); 
        init_dsol = xi_min*init_dsol;
        
        % GMRES solve
        [free_dsol,~,~,~,resvec] = gmres(free_JAC_MAT,free_RES,[],gmres_tol,maxit,L_pre,U_pre,init_dsol);

        % drop first column of Q and U and re-establish QR factorization
        if (U_d == U_dim)
            [Q, U, R] = dropFirstCol(Q,U,R);
            U_d = U_d - 1;
        end
        
        % CGS-2 low synch
        orthg_RES = free_JAC_MAT*free_dsol;
        Q(:,U_d+1) =  orthg_RES; 
        [Q, R(1:U_d+1,1:U_d+1), T] = CGS2_twoSynch(Q(:,1:U_d+1),R(1:U_d+1,1:U_d+1),T);
        
        % projection and normalization of U
        U(:,U_d+1) = free_dsol - U(:,1:U_d)*R(1:U_d,U_d+1);
        U(:,U_d+1) = U(:,U_d+1)/R(U_d+1,U_d+1);

    end

    if (U_dim > 0)
        U_d = U_d + 1; % update solution history counter
    end
    
    loc_count = loc_count+(length(resvec)-1); % update GMRES iteration counter
    
    exch_it = exch_it + 1; % update exchange iteration counter
    
    fulldsol(glb_fdof) = free_dsol; % update change in solution array
    
    fprintf('GMRES iteration = %d \n', loc_count);    
    fprintf('GMRES residual norm = %e \n', resvec(end));    
    if(resvec(end) < eps) % break if absolute residual is very low
        break;
    end

    dsol_diff = sqrt( sum((free_dsol - dsol_prev).^2) / size(glb_fdof,1) );
    fprintf('dsol diff = %e \n', dsol_diff);
    if(dsol_diff < dsol_tol) % check if dsol is still changing
        break
    end

    glb_sol_tmp(glb_fld_dof) = glb_sol_np1(glb_fld_dof) - fulldsol(glb_fld_dof);

    % determine solution at fringe points for all meshes
    glb_sol_tmp = fringe_interpolation({mesh_obj1,mesh_obj2}, {donor_map1,donor_map2}, ...
                                        glb_sol_tmp, {nd_dof_map1,nd_dof_map2}, ov_info);

    L2 = [L2, comp_L2_err(pp,glb_sol_tmp,mesh_obj1,mesh_obj2,0)];
    fprintf('\n\n');
    
end

gmres_iter_count = gmres_iter_count + loc_count;
fprintf('Total GMRES iterations = %d \n\n', gmres_iter_count);

end
 
% Drop the first column of the QR factorization, zeroing the last columns of Q
% and R.
function [Q, U, R] = dropFirstCol(Q, U, R)

dim = size(Q, 2);

R = R(:, 2:end);

for i = 1:1:(dim - 1)
    G = planerot(R(i:(i + 1), i));
    R(i:(i + 1), :) = G*R(i:(i + 1), :);
    
    Q(:, i:(i + 1)) = Q(:, i:(i + 1))*G';
    
    U(:, i:(i + 1)) = U(:, i:(i + 1))*G';
end

Q(:, dim) = 0;
U(:, dim) = 0;

R = [R zeros(dim, 1)];
end

% Classical Gram Schmidt algorithm  
function [Q, R, T] = CGS2_twoSynch(Q, R, T)

j = size(Q,2);

if (j>1)     
        tmp = Q(:, 1:j-1)'*Q(:, j-1:j);
        T(1:j-1, j-1)  = tmp(:,1);
        a = tmp(:,2);
        
        L(1:j-1, 1:j-1) = tril(T(1:j-1, 1:j-1),-1);
        
        R(1:j-1, j) = (eye(j-1,j-1) - L(:, 1:j-1) - L(:, 1:j-1)')*a;
        
        Q(:,j) = Q(:,j) - Q(:, 1:j-1)*R(1:j-1,j);
 
        Q(:,j) = Q(:,j) - Q(:,1:j-1)*(Q(:,1:j-1)'*Q(:,j));
             
end

R(j, j) = norm(Q(:,j));
Q(:, j) =  Q(:, j)/norm(Q(:,j));

end