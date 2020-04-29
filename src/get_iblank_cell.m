function iblank_cell = get_iblank_cell(mesh_def,iblank)

% get number of cells
box = mesh_def('dim');
h   = mesh_def('size');
omega_x = box(1,2) - box(1,1);
omega_y = box(2,2) - box(2,1);
cells_x = floor(omega_x/h(1));
cells_y = floor(omega_y/h(2));

% initialize iblank array
iblank_cell = ones(cells_x*cells_y,1);

for j=1:cells_y
    for i = 1:cells_x
        cell_id = (j-1)*cells_x + i;
        
        count = 0;
        for ndj=j:j+1
            for ndi=i:i+1
                nd_id = (ndj-1)*(cells_x+1) + ndi;   
                
                if(iblank(nd_id) == -1)
                    count = count+1;
                end 
            end
        end
        
        if(count == 4) % if all nodes are fringe, cell is fringe
            iblank_cell(cell_id) = -1;
        end
        
    end
end

end