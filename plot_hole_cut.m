function plot_hole_cut(mesh1, mesh2, ib1, ib2, don1, don2)

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
subplot(1,3,1)
hold on
surf(coord_x_1,coord_y_1,surf1,'edgecolor','black');
plot3(coords1(ib1==-1,1),coords1(ib1==-1,2),ib1(ib1==-1),'.','Color',[0, 0.5, 0],'markersize',25);
surf(coord_x_2,coord_y_2,surf2,'edgecolor','red');
plot3(coords2(ib2==-1,1),coords2(ib2==-1,2),ib2(ib2==-1),'.b','markersize',25);
alpha(0.0);
title('Fringe points, both meshes');
colorbar;
set(gcf,'color','w');
set(gca, 'FontSize', 18);
view(2);

subplot(1,3,2)
hold on
surf(coord_x_1,coord_y_1,surf1,'edgecolor','black');
plot3(coords1(ib1==-1,1),coords1(ib1==-1,2),ib1(ib1==-1),'.','Color',[0, 0.5, 0],'markersize',25);
surf(coord_x_2,coord_y_2,surf2,'edgecolor','red');
plot3(coords2(don2,1),coords2(don2,2),ib2(don2),'.m','markersize',25);
alpha(0.0);
title('Mesh 1 fringe with mesh 2 donor');
colorbar;
set(gcf,'color','w');
set(gca, 'FontSize', 18);
view(2);

subplot(1,3,3)
hold on
surf(coord_x_1,coord_y_1,surf1,'edgecolor','black');
plot3(coords1(don1,1),coords1(don1,2),ib1(don1),'.','Color',[0.75, 0.75, 0],'markersize',25);
surf(coord_x_2,coord_y_2,surf2,'edgecolor','red');
plot3(coords2(ib2==-1,1),coords2(ib2==-1,2),ib2(ib2==-1),'.b','markersize',25);
alpha(0.0);
title('Mesh 2 fringe with mesh 1 donor');
colorbar;
set(gcf,'color','w');
set(gca, 'FontSize', 18);
view(2);

end