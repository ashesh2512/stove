function [glb_res,glb_jac] = Build_res_jac_coupled(glb_mesh,glb_donor_map,glb_iblank, ....
                                                   glb_sol_np1,glb_sol_n,glb_sol_nm1,glb_nd_dof_map, ...
                                                   ov_info,time_info,pp)
% assemble coupled residual and jacobian correpesonding to both grids

%% pre-processing

num_grids    = ov_info('num grids'); % nmuber of grids
ndof_nd      = pp('dof per node'); % number of dofs per node
jac_srf_size = (2*ndof_nd)^2; % size of jacobain contribution from each scs

num_dofs = zeros(num_grids,1);
tot_dofs = size(glb_sol_np1,1); % total dof counter
vec_size = 0; % vector size counter
count    = 0; % counter to keep track of IVEC, JVEC, VVEC, indexing

% loop over number of grids
for ig = 1:num_grids
     
    coords = glb_mesh{ig}{2,1}; % extract coordinates
    
    edge_topo_x = glb_mesh{ig}{2,2}; % extract edge topologies
    
    edge_topo_y = glb_mesh{ig}{2,4}; % extract edge areas
    
    num_dofs(ig) = ndof_nd*size(coords,1); % number of dofs on this mesh
    
    donor_map = glb_donor_map{ig}; % donor map array for current grid
    
    % determine size for compressed sparse column format vector
    vec_size = vec_size ...
             + jac_srf_size*(size(edge_topo_x,1) + size(edge_topo_y,1)) ...
             + ndof_nd*(size(coords,1) - size(donor_map,1)) ...
             + ndof_nd*(sum(cellfun('length',donor_map(:,2))) + size(donor_map,1));
    
end

% initialize global residual array
glb_res = zeros(tot_dofs,1);

% compressed sparse column format vector components
IVEC = zeros(vec_size,1);
JVEC = zeros(vec_size,1);
VVEC = zeros(vec_size,1);

%% loop over all grids
for ig = 1:num_grids
    
    % extract relevant data for this mesh grid
    coords = glb_mesh{ig}{2,1}; % extract coordinates 

    edge_topo_x = glb_mesh{ig}{2,2}; % x-edge topology
    edge_area_x = glb_mesh{ig}{2,3}; % x-edge area

    edge_topo_y = glb_mesh{ig}{2,4}; % y-edge topology
    edge_area_y = glb_mesh{ig}{2,5}; % y-edge area

    nd_vol = glb_mesh{ig}{2,6}; % extract nodal volumes
    
    nd_dof_map = glb_nd_dof_map{ig}; % extract node to dof map
    
    iblank = glb_iblank{ig}; % iblank array for current grid

    donor_map = glb_donor_map{ig}; % donor map array for fringe nodes

    % extract donor mesh info
    donor_grid = ov_info(strcat('mesh',num2str(ig),' donor'));
    
    donor_coords = glb_mesh{donor_grid}{2,1}; % extract donor mesh coordinates
    
    donor_nd_dof_map = glb_nd_dof_map{donor_grid}; % extract node to dof map
    
    %% loop over all edges aligned with x for surface contributions
    for ex = 1:size(edge_topo_x,1)
        
        edge_nd = edge_topo_x(ex,:); % extract nodes belonging to the edge
        
        % extract only x-coordinates for the current edge
        edge_coord = coords(edge_nd,1)';
        
        % define normal to the edge in +ve x direction
        edge_nrml = [0; 1];
        
        % extract dof ids associated with nodes of current edge
        edge_dofs = reshape(nd_dof_map(edge_nd,:)',[],1);
        
        % extract field dofs
        field_dofs = reshape(nd_dof_map(edge_nd(iblank(edge_nd)==1),:)',[],1);

        % compute residual and jacobian contribution from surface
        % contributions
        [res_srf, jac_srf] = compute_res_jac_srf( glb_sol_np1(edge_dofs), edge_coord, edge_area_x(ex), edge_nrml, pp );
        
        % sum into global reidual and jacobian for field dofs
        glb_res(field_dofs) = glb_res(field_dofs) + res_srf(ismember(edge_dofs,field_dofs));
        
        % store jac_srf in CSC format
        for r = 1:length(edge_dofs)
            for c = 1:length(edge_dofs)
                count = count + 1;
                IVEC(count,1) = edge_dofs(r);
                JVEC(count,1) = edge_dofs(c);
                
                % add jacobian contribution only for field points
                if( ismember(edge_dofs(r),field_dofs) ) 
                    VVEC(count,1) = jac_srf(r,c);
                end
            end
        end
        
    end
    
    %% loop over all edges aligned with y for surface contributions
    for ey = 1:size(edge_topo_y,1)
        
        edge_nd = edge_topo_y(ey,:); % extract nodes belonging to the edge
        
        % extract only y-coordinates for the current edge
        edge_coord = coords(edge_nd,2)';
        
        % define normal to the edge in +ve x direction
        edge_nrml = [0; 1];
        
        % extract dof ids associated with nodes of current edge
        edge_dofs = reshape(nd_dof_map(edge_nd,:)',[],1);

        % extract field dofs
        field_dofs = reshape(nd_dof_map(edge_nd(iblank(edge_nd)==1),:)',[],1);
        
        % compute residual and jacobian contribution from surface
        % contributions
        [res_srf, jac_srf] = compute_res_jac_srf( glb_sol_np1(edge_dofs), edge_coord, edge_area_y(ey), edge_nrml, pp );
        
        % sum into global reidual and jacobian for field dofs
        glb_res(field_dofs) = glb_res(field_dofs) + res_srf(ismember(edge_dofs,field_dofs));

        % store jac_srf in CSC format
        for r = 1:length(edge_dofs)
            for c = 1:length(edge_dofs)
                count = count + 1;
                IVEC(count,1) = edge_dofs(r);
                JVEC(count,1) = edge_dofs(c);
                
                % add jacobian contribution only for field points
                if( ismember(edge_dofs(r),field_dofs) ) 
                    VVEC(count,1) = jac_srf(r,c);
                end
            end
        end
        
    end
    
    %% loop over all nodes for volume contribution
    for nd = 1:size(coords,1)
        
        % skip loop for fringe nodes
        if( iblank(nd) == -1 )
            continue;
        end
        
        % extract only coordinates for the current node
        coord = coords(nd,:);
        
        % extract dof ids associated with current node
        nd_dofs = nd_dof_map(nd,:)';
   
        % compute residual and jacobian contribution from volumetric
        % contributions
        [res_vol, jac_vol] = compute_res_jac_vol(glb_sol_np1(nd_dofs), glb_sol_n(nd_dofs), glb_sol_nm1(nd_dofs), ...
                                                 coord, nd_vol(nd), time_info, pp);
        
        % sum into global reidual and jacobian
        glb_res(nd_dofs) = glb_res(nd_dofs) + res_vol;

        % store jac_vol in CSC format
        for r = 1:length(nd_dofs)
            for c = 1:length(nd_dofs)
                count = count + 1;
                IVEC(count,1) = nd_dofs(r);
                JVEC(count,1) = nd_dofs(c);
                VVEC(count,1) = jac_vol(r,c);
            end
        end
        
    end
    
    %% loop over all fringe nodes for overset contribution
    for ifr = 1:size(donor_map,1)
        
        frng_nd   = donor_map{ifr,1}; % fringe node
        donor_nds = donor_map{ifr,2}; % extract donor node ids
        
        frng_coords     = coords(frng_nd,:); % extract fringe node coordinates
        donor_nd_coords = donor_coords(donor_nds,:); % extract donor node coordinates
        
        frng_nd_dofs  = reshape(nd_dof_map(frng_nd,:)',[],1); % extract fringe node dofs
        donor_nd_dofs = reshape(donor_nd_dof_map(donor_nds,:)',[],1); % extract donor node dofs 
        
        % compute residual and jacobian contribution from interpolation
        % contributions
        coeff = compute_frg_coeff(frng_coords,donor_nd_coords,ov_info);
                
        % sum into global reidual and jacobian
        for r = 1:length(frng_nd_dofs)
            % sum into residual = frng_sol - N*don_sol
            glb_res(frng_nd_dofs(r)) = glb_res(frng_nd_dofs(r)) ... 
                                     + glb_sol_np1(frng_nd_dofs(r)) - coeff'*glb_sol_np1(donor_nd_dofs);

            % jacobian contribution from d(frng_sol)/d(frng_sol)
            count = count + 1;
            IVEC(count,1) = frng_nd_dofs(r);
            JVEC(count,1) = frng_nd_dofs(r);
            VVEC(count,1) = 1.0;
            
            % jacobian contribution from d(-N*don_sol)/d(don_sol)
            for c = 1:length(donor_nd_dofs)
                count = count + 1;
                IVEC(count,1) = frng_nd_dofs(r);
                JVEC(count,1) = donor_nd_dofs(c);
                VVEC(count,1) = -coeff(c);
            end
        end

    end
    
end

% create sparse matrix out of IVEC,JVEC,VVEC
glb_jac = full(sparse(IVEC,JVEC,VVEC,tot_dofs,tot_dofs));

end