function [nb_ib, bg_ib] = cut_hole(nb_mesh, bg_mesh, nb_bc, ov_info)
% perform holecutting based on mandatory fringe points on near body mesh
% create iblank arrays with field points marked as 1 and fringe points
% marked as -1

% extract edge topologies
nb_edge_topo_x = nb_mesh{2,2};
nb_edge_topo_y = nb_mesh{2,4};

% extract sidesets
bottom_ss = nb_mesh{2,7};
right_ss  = nb_mesh{2,8};
top_ss    = nb_mesh{2,9};
left_ss   = nb_mesh{2,10};

% extract mesh coordinates 
nb_coords = nb_mesh{2,1};
bg_coords = bg_mesh{2,1};

% declare iblank arrays
nb_ib = ones(size(nb_coords,1),1);
bg_ib = ones(size(bg_coords,1),1);

% extract number of mandatory fringe points along each overset boundary
% constant for entire mesh
mndtry_frng = ov_info('mandatory frng');

% extract minimum desired overset boundary-speicific overlap arranged 
% counter-clockwise
overlap = ov_info('overlap');

tol = 1e-14; % tolerance for comparison of coordinates

%% collect all fringe points correponding to bottom overset boundary

if (nb_bc('bottom') == "overset") % assuming mesh was constructed bottom to top
    frng_nds_nb = reshape(nb_edge_topo_x(1 : mndtry_frng*length(bottom_ss), :), [], 1);
    nb_ib(frng_nds_nb) = -1;
    
    % blankly mark everything that is to above fringe + overlap as fringe
    % for background
    frng_nds_bg = bg_coords(:,2) >= (max(nb_coords(frng_nds_nb,2))+overlap(1)-tol);
    bg_ib(frng_nds_bg) = -1;
end

%% collect all fringe points correponding to right overset boundary

if (nb_bc('right') == "overset") % assuming mesh was constructed left to right
    frng_nds_nb = reshape(nb_edge_topo_y(right_ss(end)-mndtry_frng*length(right_ss)+1 : right_ss(end), :), [], 1);
    nb_ib(frng_nds_nb) = -1;
    
    bg_ib_tmp = bg_ib; % store current iblank configuration
    
    % blankly mark everything that is to the left of fringe - overlap as
    % fringe for background
    frng_nds_bg = bg_coords(:,1) <= (min(nb_coords(frng_nds_nb,1))-overlap(2)+tol);
    bg_ib(frng_nds_bg) = -1;
        
    bg_ib(bg_ib_tmp == 1) = 1; % restore field points

    % blankly mark everything that isn't frng_nds_bg as field points
    field_nds_bg = bg_coords(:,1) >= (min(nb_coords(frng_nds_nb,1))-overlap(2)+tol);
    bg_ib(field_nds_bg) = 1;
end

%% collect all fringe points correponding to top overset boundary

if (nb_bc('top') == "overset") % assuming mesh was constructed bottom to top
    frng_nds_nb = reshape(nb_edge_topo_x(top_ss(end)-mndtry_frng*length(top_ss)+1 : top_ss(end), :), [], 1);
    nb_ib(frng_nds_nb) = -1;
    
    bg_ib_tmp = bg_ib; % store current iblank configuration
    
    % blankly mark everything that is to the bottom of fringe - overlap as
    % fringe for background
    frng_nds_bg = bg_coords(:,2) <= (min(nb_coords(frng_nds_nb,2))-overlap(3)+tol);
    bg_ib(frng_nds_bg) = -1;
        
    bg_ib(bg_ib_tmp == 1) = 1; % restore field points

    % blankly mark everything that isn't frng_nds_bg as field points
    field_nds_bg = bg_coords(:,2) >= (min(nb_coords(frng_nds_nb,2))-overlap(3)+tol);
    bg_ib(field_nds_bg) = 1;
end

%% collect all fringe points correponding to left overset boundary

if (nb_bc('left') == "overset") % assuming mesh was constructed left to right
    frng_nds_nb = reshape(nb_edge_topo_y(1 : mndtry_frng*length(left_ss),:),[],1);
    nb_ib(frng_nds_nb) = -1;
    
    bg_ib_tmp = bg_ib; % store current iblank configuration
    
    % blankly mark everything that is to the left of fringe - overlap as
    % fringe for background
    frng_nds_bg = bg_coords(:,1) >= (max(nb_coords(frng_nds_nb,1))+overlap(4)-tol);
    bg_ib(frng_nds_bg) = -1;
        
    bg_ib(bg_ib_tmp == 1) = 1; % restore field points

    % blankly mark everything that isn't frng_nds_bg as field points
    field_nds_bg = bg_coords(:,1) <= (max(nb_coords(frng_nds_nb,1))+overlap(4)-tol);
    bg_ib(field_nds_bg) = 1;  
end

end