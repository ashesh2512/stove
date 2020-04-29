function mesh_obj_cell = meshgen_cell(box, h)
% structured mesh generator for cell-centered finite volume. Ordering
% follows a left to right, bottom to top pattern.
 
xlim = box(1,:);
ylim = box(2,:);

% 1D mesh along x
omega_x  = xlim(2) - xlim(1);
ncells_x = floor(omega_x/h(1));
coord_x  = linspace(xlim(1)+h(1)/2, xlim(2)-h(1)/2, ncells_x)';

% 1D mesh along y
omega_y  = ylim(2) - ylim(1);
ncells_y = floor(omega_y/h(2));
coord_y  = linspace(ylim(1)+h(2)/2, ylim(2)-h(2)/2, ncells_y)';

% intialize arrays
coord_list  = zeros(ncells_x*ncells_y,2);

%% loop along y-axis and construct edges left to right
for iy = 1:ncells_y
    
    % fill coordinate list
    cell_ids_x = (((iy-1)*ncells_x+1):(iy*ncells_x))';
    coord_list(cell_ids_x,:) = [coord_x coord_y(iy)*ones(ncells_x,1)];
end

%% store variables of interest
% any change in ordering here should reflect access throughout the code

mesh_obj_cell{1, 1} =  "coord_list"; mesh_obj_cell{2, 1} = coord_list;

end