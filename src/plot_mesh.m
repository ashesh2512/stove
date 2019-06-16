function plot_mesh(mesh1, mesh2)

% 1D mesh along x
coords1   = mesh1{2,1};
coord_x_1 = unique(coords1(:,1));
coord_y_1 = unique(coords1(:,2));
[X1,Y1]   = meshgrid(coord_x_1,coord_y_1);
surf1     = 0.*X1 + 0.*Y1;

% 1D mesh along y
coords2   = mesh2{2,1};
coord_x_2 = unique(coords2(:,1));
coord_y_2 = unique(coords2(:,2));
[X2,Y2]   = meshgrid(coord_x_2,coord_y_2);
surf2     = 0.*X2 + 0.*Y2;

figure()
hold on
surf(coord_x_1,coord_y_1,surf1,'edgecolor','black','facecolor','white');
surf(coord_x_2,coord_y_2,surf2,'edgecolor','red','facecolor','white');
alpha(0.0);
title('overset mesh');
colorbar;
set(gcf,'color','w');
set(gca, 'FontSize', 18);
view(2);

end