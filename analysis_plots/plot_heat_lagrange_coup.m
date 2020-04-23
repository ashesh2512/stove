clc; close all; clear all;

%% problem properties
% pp = containers.Map({            'prblm', 'dof per node', 'conductivity'}, ...
%                      { "steady heat MMS",              1,              1});
%                  
% % background mesh
% box1 = [-2,2; -2,2];
% h1   = [ 0.2,  0.2]/2/2/2;
% 
% % boundary condition map - bottom, right, top, left
% bc1 = containers.Map({   'bottom',     'right',       'top',      'left'}, ...
%                      {"dirichlet", "dirichlet", "dirichlet", "dirichlet"});
% % near body mesh
% box2 = [-1 1; -1 1];
% h2   = [ 0.2,  0.2]/2/2/2;
% 
% % boundary condition map - bottom, right, top, left
% bc2 = containers.Map({ 'bottom',   'right',     'top',    'left'}, ...
%                      {"overset", "overset", "overset", "overset"});

% ov_info = containers.Map({ 'num grids', 'mesh1 donor', 'mesh2 donor', 'mandatory frng', 'overlap', 'donor_grid', 'intrp_type', 'intrp_order'}, ...
%                          {  2, 2, 1, 2, [4*0.1, 4*0.1, 4*0.1, 4*0.1], "tensor", "lagrange", 1 });
% 
% % convergence parameters
% N_iters    =    10; % Maximum number of Newton steps
% resnrmdrop = 1e-09; % Newton convergence criteria

del    = 0.2;
h      = [del del/2 del/4 del/8];
L2     = [4.598126e-02 1.113857e-02 2.782573e-03 6.979062e-04];
y_vals = 15*h.^2;


figure(1)
clf
hold on
plot(h,y_vals,  '-k','LineWidth', 2);
plot(h,L2,'--*r','LineWidth', 2,'MarkerSize',20);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('$h$','Interpreter','latex'); ylabel('$L_2$','Interpreter','latex');
% title('coincident meshes','Interpreter','latex');
legend('slope=2', 'single mesh', 'Interpreter','latex');
set(gcf,'color','w');
set(gca, 'FontSize', 24);

%% problem properties               
% % background mesh
% box1 = [-2,2; -2,2];
% h1   = [ 0.2,  0.2];
%
% near body mesh
% box2 = [-1.13625,0.86375; -1.13625,0.86375];
% h2   = [             0.2,              0.2];
%
% ov_info = containers.Map({ 'num grids', 'mesh1 donor', 'mesh2 donor', 'mandatory frng', 'overlap', 'donor grid', 'intrp order'}, ...
%                          {  2, 2, 1, 2, [3*h2(1), 3*h2(1), 3*h2(1), 3*h2(1)], "tensor", 1 });

del    = 0.2;
h      = [del del/2 del/4 del/8];
L2_p1  = [6.526605e-02 9.348554e-03 6.497774e-03 5.626061e-03];
L2_p2  = [4.684465e-02 8.873482e-03 3.253635e-03 6.685104e-04];
L2_p3  = [4.962786e-02 1.088594e-02 2.814165e-03 7.110412e-04];
y_vals1= 1*h;
y_vals2= 0.5*h.^2;

figure(2)
subplot(1,2,1)
clf
hold on
plot(h,y_vals1,'-m','LineWidth', 2);
plot(h,y_vals2,'-k','LineWidth', 2);
plot(h,L2_p1,'--sr','LineWidth', 2,'MarkerSize',20);
plot(h,L2_p2,'--^g','LineWidth', 2,'MarkerSize',20);
plot(h,L2_p3,'--ob','LineWidth', 2,'MarkerSize',20);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('$h$','Interpreter','latex'); ylabel('$L_2$','Interpreter','latex');
legend('slope=1', 'slope=2', 'linear Lagrange', 'quadratic Lagrange', 'cubic Lagrange', 'Interpreter','latex');
set(gcf,'color','w');
set(gca, 'FontSize', 24);

gap = [7.362500e-01, 6.637500e-01, 6.637500e-01, 7.362500e-01;
       3.362500e-01, 3.637500e-01, 3.637500e-01, 3.362500e-01;
       1.862500e-01, 1.637500e-01, 1.637500e-01, 1.862500e-01;
       8.625000e-02, 8.875000e-02, 8.875000e-02, 8.625000e-02];

figure(3)
clf
subplot(2,2,1)
plot(h, gap(:,1),'--*r','LineWidth', 2)
xlabel('h','Interpreter','latex'); ylabel('fringe gap','Interpreter','latex');
title('bottom boundary','Interpreter','latex');
set(gcf,'color','w');
set(gca, 'FontSize', 16);

subplot(2,2,2)
plot(h, gap(:,2),'--*r','LineWidth', 2)
xlabel('h','Interpreter','latex'); ylabel('fringe gap','Interpreter','latex');
title('right boundary','Interpreter','latex');
set(gcf,'color','w');
set(gca, 'FontSize', 16);

subplot(2,2,3)
plot(h, gap(:,3),'--*r','LineWidth', 2)
xlabel('h','Interpreter','latex'); ylabel('fringe gap','Interpreter','latex');
title('top boundary','Interpreter','latex');
set(gcf,'color','w');
set(gca, 'FontSize', 16);

subplot(2,2,4)
plot(h, gap(:,4),'--*r','LineWidth', 2)
xlabel('h','Interpreter','latex'); ylabel('fringe gap','Interpreter','latex');
title('left boundary','Interpreter','latex');
set(gcf,'color','w');
set(gca, 'FontSize', 16);

fig = figure(2);
set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig,'temp_name','-dpdf','-r0')

%% problem properties               
% near body mesh 2 fringe gap: bottom: 7.362500e-01, right: 6.637500e-01, top: 6.637500e-01, left:7.362500e-01 

% ov_info = containers.Map({ 'num grids', 'mesh1 donor', 'mesh2 donor', 'mandatory frng', 'overlap', 'donor grid', 'intrp order'}, ...
%                          {  2, 2, 1, 2, [3*h2(1), 3*h2(1), 3*h2(1), 3*h2(1)], "tensor", 1 });

% ov_info = containers.Map({ 'num grids', 'mesh1 donor', 'mesh2 donor', 'mandatory frng', 'overlap', 'donor grid', 'intrp order'}, ...
%                          {  2, 2, 1, 3, [7*h2(1), 6*h2(1), 6*h2(1), 7*h2(1)], "tensor", 1 });

% ov_info = containers.Map({ 'num grids', 'mesh1 donor', 'mesh2 donor', 'mandatory frng', 'overlap', 'donor grid', 'intrp order'}, ...
%                          {  2, 2, 1, 5, [14*h2(1), 13*h2(1), 13*h2(1), 14*h2(1)], "tensor", 1 });

% ov_info = containers.Map({ 'num grids', 'mesh1 donor', 'mesh2 donor', 'mandatory frng', 'overlap', 'donor grid', 'intrp order'}, ...
%                          {  2, 2, 1, 9, [29*h2(1), 26*h2(1), 26*h2(1), 29*h2(1)], "tensor", 1 });
                     
del    = 0.2;
h      = [del del/2 del/4 del/8];
L2_p1  = [6.526605e-02 1.619219e-02 3.762973e-03 1.023410e-03];
y_vals1= 5*h;
y_vals2= 0.1*h.^2;

figure(4)
clf
hold on
plot(h,y_vals2,'-k','LineWidth', 2);
plot(h,L2_p1,'--*r','LineWidth', 2,'MarkerSize',20);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('$h$','FontSize',18,'Interpreter','latex'); ylabel('$L_2$','FontSize',18,'Interpreter','latex');
legend('slope=2', 'linear Lagrange', 'Interpreter','latex');
% title('constant overlap','Interpreter','latex');
set(gcf,'color','w');
set(gca, 'FontSize', 18);
