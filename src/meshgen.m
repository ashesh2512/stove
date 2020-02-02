function mesh_obj = meshgen(box, h)
% structured mesh generator for vertex-centered finite volume. Ordering
% follows a left to right, bottom to top pattern.
%
% Input:  box  - coordinates of mesh extrema
%         h    - element size of mesh
%         bc   - map containing bc for each boundary
%
% Output: mesh_obj - mesh object
 
xlim = box(1,:);
ylim = box(2,:);

% 1D mesh along x
omega_x  = xlim(2) - xlim(1);
nnodes_x = floor(omega_x/h(1)) + 1;
coord_x  = linspace(xlim(1), xlim(2), nnodes_x)';

% 1D mesh along y
omega_y  = ylim(2) - ylim(1);
nnodes_y = floor(omega_y/h(2)) + 1;
coord_y  = linspace(ylim(1), ylim(2), nnodes_y)';

% intialize arrays
coord_list  = zeros(nnodes_x*nnodes_y,2);

edge_topo_x = zeros((nnodes_x-1)*(nnodes_y),2);
edge_area_x = ones((nnodes_x-1)*nnodes_y,1)*h(2); % area of scs normal to the edge

edge_topo_y = zeros((nnodes_y-1)*(nnodes_x),2);
edge_area_y = ones((nnodes_y-1)*nnodes_x,1)*h(1); % area of scs normal to the edge

nd_vol = ones(nnodes_x*nnodes_y,1)*0.0625; % volume of scv
                                           % multiplication factor of
                                           % 0.0625 is more of a coeff
                                           % that is doubled when an edge
                                           % shares the same node in a
                                           % Quad4 element

%% loop along y-axis and construct edges left to right
for iy = 1:nnodes_y
    
    % fill coordinate list
    node_ids_x = (((iy-1)*nnodes_x+1):(iy*nnodes_x))';
    coord_list(node_ids_x,:) = [coord_x coord_y(iy)*ones(nnodes_x,1)];

    % create edge topology for x-direction
    edge_ids_x = (((iy-1)*(nnodes_x-1)+1):(iy*(nnodes_x-1)))';
    edge_topo_x(edge_ids_x,:) = [node_ids_x(1:end-1) node_ids_x(2:end)];
end

%% loop along x-axis and construct edges bottom to top
for ix = 1:nnodes_x

    % create edge topology for y-direction
    node_ids_y = (ix:nnodes_x:(nnodes_x*nnodes_y))';
    edge_ids_y = (((ix-1)*(nnodes_y-1)+1):(ix*(nnodes_y-1)))';
    edge_topo_y(edge_ids_y,:) = [node_ids_y(1:end-1) node_ids_y(2:end)];
end

%% compute the surface and volume areas 
tol = 1e-14; % tolerance for comparison of coordinates

% loop over all edges along x-axis left to right
for i = 1:size(edge_topo_x,1)
    
    nd = edge_topo_x(i,:);
    
    % edges along boundaries have scs with half the area
    if( norm(coord_list(nd,2)-ylim(1)) < tol || norm(coord_list(nd,2)-ylim(2)) < tol )
        edge_area_x(i) = edge_area_x(i)/2;
    end
    
    nd_vol(nd) = nd_vol(nd)*2; % volume coefficient update for every edge share
end

% loop over all edges along y-axis bottom to top
for i = 1:size(edge_topo_y,1)
    
    nd = edge_topo_y(i,:);
    
    % edges along boundaries have scs with half the area
    if( norm(coord_list(nd,1)-xlim(1)) < tol || norm(coord_list(nd,1)-xlim(2)) < tol )
        edge_area_y(i) = edge_area_y(i)/2;
    end
    
    nd_vol(nd) = nd_vol(nd)*2; % volume coefficient update for every edge share
end

nd_vol = nd_vol.*(h(1)*h(2)); % scale volume with coefficients computed above

%% store side sets in terms of edge ids
% bottom and top sidesets stored in terms of edge_topo_x
% left and right sidesets stored in terms of edge_topo_y

bottom_ss = (1:(nnodes_x-1))';
right_ss  = (((nnodes_y-1)*(nnodes_x-1)+1):(nnodes_y-1)*(nnodes_x))';
top_ss    = (((nnodes_x-1)*(nnodes_y-1)+1):(nnodes_x-1)*(nnodes_y))';
left_ss   = (1:(nnodes_y-1))';

%% store variables of interest
% any change in ordering here should reflect access throughout the code

mesh_obj{1, 1} =  "coord_list"; mesh_obj{2, 1} = coord_list;
mesh_obj{1, 2} = "edge_topo_x"; mesh_obj{2, 2} = edge_topo_x;
mesh_obj{1, 3} = "edge_area_x"; mesh_obj{2, 3} = edge_area_x;
mesh_obj{1, 4} = "edge_topo_y"; mesh_obj{2, 4} = edge_topo_y;
mesh_obj{1, 5} = "edge_area_y"; mesh_obj{2, 5} = edge_area_y;
mesh_obj{1, 6} =      "nd_vol"; mesh_obj{2, 6} = nd_vol;
mesh_obj{1, 7} =   "bottom_ss"; mesh_obj{2, 7} = bottom_ss;
mesh_obj{1, 8} =    "right_ss"; mesh_obj{2, 8} = right_ss;
mesh_obj{1, 9} =      "top_ss"; mesh_obj{2, 9} = top_ss;
mesh_obj{1,10} =     "left_ss"; mesh_obj{2,10} = left_ss;

end