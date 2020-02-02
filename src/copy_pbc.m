function glb_sol = copy_pbc(num_grid,pnd,nd_dof_map,glb_sol)

for ig = 1:num_grid
    
    pnd_ig = pnd{ig}; % exctract periodic node pairs for current grid
    if(isempty(pnd_ig))
        continue
    end
    
    nd_dof_map_ig = nd_dof_map{ig}; % exctract node to dof maps for current grid
    
    for nd = 1:size(pnd_ig,1)
        
        mas_nd_dofs = nd_dof_map_ig(pnd_ig(nd,1),:)'; % master node dofs
        slv_nd_dofs = nd_dof_map_ig(pnd_ig(nd,2),:)'; % slave node dofs
        
        % copy solution from master dofs to slave dofs
        for idof = 1:length(mas_nd_dofs)
            glb_sol(slv_nd_dofs(idof)) = glb_sol(mas_nd_dofs(idof));
        end
    
    end
    
end

end