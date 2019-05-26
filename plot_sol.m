function plot_sol(mesh1,mesh2,sol,nd_dof_map1,nd_dof_map2)

for idof = 1:size(nd_dof_map1,2)
    % mesh 1 solution
    coords1   = mesh1{2,1};
    coord_x_1 = unique(coords1(:,1));
    coord_y_1 = unique(coords1(:,2));
    surf1     = reshape(sol(nd_dof_map1(:,idof)),[length(coord_x_1),length(coord_y_1)])';
    
    % mesh 2 solution
    coords2   = mesh2{2,1};
    coord_x_2 = unique(coords2(:,1));
    coord_y_2 = unique(coords2(:,2));
    surf2     = reshape(sol(nd_dof_map2(:,idof)),[length(coord_x_2),length(coord_y_2)])';
    
    figure()
    hold on
    surf(coord_x_1,coord_y_1,surf1,'edgecolor','none');
    surf(coord_x_2,coord_y_2,surf2,'edgecolor','none');
    title('overset solution');
    colorbar;
    set(gcf,'color','w');
    set(gca, 'FontSize', 18);
    view(2); 
end

end