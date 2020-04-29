function sol_out = shift_solution(mesh,ib,sol_in_nd,sol_in_cell,nd_dof_map,cell_dof_map,averaging,ov_info)

num_grids = length(mesh); % nmuber of grids
ndof = size(nd_dof_map{1},2); % number of dofs per node/cell

for ig = 1:num_grids
    
    if(ov_info('background grid')~=ig)
        continue
    end
    
    % get number of cells
    coords = mesh{ig}{2,1};
    cells_x = size(unique(coords(:,1)),1)-1;
    cells_y = size(unique(coords(:,2)),1)-1;
    
    switch averaging
        case "cell" % node to cell averaging
            sol_out = sol_in_cell;
            
            for j=1:cells_y
                for i = 1:cells_x
                    cell_id = (j-1)*cells_x + i;
                    
                    sum = zeros(ndof,1);
                    for ndj=j:j+1
                        for ndi=i:i+1
                            nd_id = (ndj-1)*(cells_x+1) + ndi;
                            sum = sum + sol_in_nd(nd_dof_map{ig}(nd_id,:));
                        end
                    end
                    sol_out(cell_dof_map{ig}(cell_id,:)) = sum./4;
                end
            end
            
        case "node" % cell to node averaging
            sol_out = sol_in_nd;
            
            for ndj=1:(cells_y+1)
                for ndi = 1:(cells_x+1)
                    nd_id = (ndj-1)*(cells_x+1) + ndi;
                    
                    if(ib{ig}(nd_id) == -1)
                        sol_out(nd_dof_map{ig}(nd_id,:)) = 0; % reset solution at fringe
                        for j=(ndj-1):ndj
                            for i=(ndi-1):ndi
                                cell_id = (j-1)*cells_x + i;
                                sol_out(nd_dof_map{ig}(nd_id,:)) = sol_out(nd_dof_map{ig}(nd_id,:)) ...
                                                                 + sol_in_cell(cell_dof_map{ig}(cell_id,:))./4;
                            end
                        end
                    end
                    
                end
            end
            
        otherwise
            error('Invalid averaging type');
    end
    
end

end