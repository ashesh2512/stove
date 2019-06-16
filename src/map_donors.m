function [donor_map,bg_don_nds] = map_donors(nb_mesh, bg_mesh, nb_ib, bg_ib, ov_info)
% map donors in background mesh to  holecutting based on mandatory fringe 
% points on near body mesh create iblank arrays with field points marked as
% 1 and fringe points marked as -1

% determine donor grid
switch ov_info('donor grid')
    case "tensor"
        [donor_map, bg_don_nds] = get_nbrs_tensor(nb_mesh{2,1}, bg_mesh{2,1}, nb_ib, ov_info('intrp order'));
    
    case "radial"
        [donor_map, bg_don_nds] = get_nbrs_radial(nb_mesh{2,1}, bg_mesh{2,1}, nb_ib, ov_info('intrp radius'));

    otherwise
        error('Do not recognize donor grid type; check overset info in driver');
end

% make sure the fringe nodes from each mesh aren't donating to any fringe
% nodes on the other mesh
if(~isempty( intersect(find(bg_ib==-1), bg_don_nds) ))
    error('Fatal Error! Found fringe nodes from backgound mesh donating to near body mesh');
end

end

%% get neighbors arranged on a tensor grid

function [donor_map,bg_don_nds] = get_nbrs_tensor(nb_coords, bg_coords, nb_ib, p)

bg_xcoords = unique(bg_coords(:,1)); % x coordinates
bg_ycoords = unique(bg_coords(:,2)); % y coordinates

nb_frng_id = find(nb_ib == -1); % determine ids of fringe node
donor_map  = cell(length(nb_frng_id),2); % initialize donor map array

bg_don_nds = []; %initialize background donor array

for ifr=1:length(nb_frng_id)
    
    donor_map{ifr,1} = nb_frng_id(ifr); % store id of fringe node
    
     % compute distance of fringe node from nodes of overlapping mesh
    [~, sorted_id_x] = sort( abs(nb_coords(nb_frng_id(ifr),1) - bg_xcoords) );
    [~, sorted_id_y] = sort( abs(nb_coords(nb_frng_id(ifr),2) - bg_ycoords) );
    
    % determine structured donor grid (p+1xp+1,1) arranged row wise
    intrp_grid = (sort(sorted_id_y(1:p+1))-1)*length(bg_xcoords) + sort(sorted_id_x(1:p+1))';
    intrp_grid = reshape(intrp_grid',[],1);
    
    donor_map{ifr,2} = intrp_grid;
    
    % mark all donor nodes as mandatory field nodes
    bg_don_nds = [bg_don_nds; donor_map{ifr,2}];
end

bg_don_nds = unique(bg_don_nds); % remove duplicate entries

end

%% get neighbors arranged on a radial grid

function [donor_map,bg_don_nds] = get_nbrs_radial(nb_coords, bg_coords, nb_ib, r)

bg_xcoords = unique(bg_coords(:,1)); % x coordinates
bg_ycoords = unique(bg_coords(:,2)); % y coordinates

% compute element length in x and y direction
bg_xh = (max(bg_xcoords) - min(bg_xcoords)) / (length(bg_xcoords)-1);
bg_yh = (max(bg_ycoords) - min(bg_ycoords)) / (length(bg_ycoords)-1);

% determine number of nodes in each direction to form search box
xnbrs = 2*ceil(r/bg_xh);
ynbrs = 2*ceil(r/bg_yh);

nb_frng_id = find(nb_ib == -1); % determine ids of fringe node
donor_map  = cell(length(nb_frng_id),2); % initialize donor map array

bg_don_nds = []; %initialize background donor array

for ifr=1:length(nb_frng_id)
    
    donor_map{ifr,1} = nb_frng_id(ifr); % store id of fringe node
    
    % compute distance of fringe node from nodes of overlapping mesh
    [~, sorted_id_x] = sort( abs(nb_coords(nb_frng_id(ifr),1) - bg_xcoords) );
    [~, sorted_id_y] = sort( abs(nb_coords(nb_frng_id(ifr),2) - bg_ycoords) );
    
    % build structured search grid based on xnbrs and ynbrs
    search_grid = (sorted_id_y(1:ynbrs)-1)*length(bg_xcoords) + sorted_id_x(1:xnbrs)';

    % determine point cloud based on search grid    
    donor_map{ifr,2} = search_grid(vecnorm(nb_coords(nb_frng_id(ifr),:) - bg_coords(search_grid,:),2,2) <= r);
    
    % mark all donor nodes as mandatory field nodes
    bg_don_nds = [bg_don_nds; donor_map{ifr,2}];
    
end

bg_don_nds = unique(bg_don_nds); % remove duplicate entries

end