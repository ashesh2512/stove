function plot_sol_cell(mesh,sol,cell_dof_map,time)

% reshape cell iblank
coords = mesh{2,1};
coord_x = unique(coords(:,1));
coord_y = unique(coords(:,2));
xh = abs(coord_x(2) - coord_x(1)); yh = abs(coord_y(2) - coord_y(1));
x_ibc = [min(coord_x)+0.5*xh,max(coord_x)-0.5*xh];
y_ibc = [min(coord_y)+0.5*yh,max(coord_y)-0.5*yh];

for idof = 1:size(cell_dof_map,2)

    solplot = reshape(sol(cell_dof_map(:,idof)),size(coord_x,1)-1,size(coord_y,1)-1);

    figure()
    image(x_ibc,y_ibc,solplot,'CDataMapping','scaled');
    colormap(parula);
    colorbar;
    title(sprintf('overset solution on background mesh at %f',time));
    set(gcf,'color','w');
    set(gca, 'FontSize', 28);
    set(gca,'YDir','normal');

end

end