clc; close all; clear all;

%% problem properties
% % problem properties
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
% box2 = [-1.13625,0.86375; -1.13625,0.86375];
% h2   = [             0.2,              0.2]/2/2/2;
% 
% % boundary condition map - bottom, right, top, left
% bc2 = containers.Map({ 'bottom',   'right',     'top',    'left'}, ...
%                      {"overset", "overset", "overset", "overset"});
%
% ov_info = containers.Map({ 'num grids', 'mesh1 donor', 'mesh2 donor', 'mandatory frng', 'overlap', 'donor grid', ...
%                            'intrp radius', 'intrp type', 'intrp shape', 'shape param', 'poly order'}, ...
%                          { 2, 2, 1, 2, [3*h2(1), 3*h2(1), 3*h2(1), 3*h2(1)], "radial", ...
%                            2.5*max(h2), "rbf", "gaussian", 1.0, 3 });
% 
% % convergence parameters
% N_iters    =    10; % Maximum number of Newton steps
% resnrmdrop = 1e-09; % Newton convergence criteria

del    = 0.2;
h      = [del del/2 del/4 del/8];
L2_pm1 = [5.473572e-02 1.031028e-02 3.263450e-03 1.435484e-03];
L2_p1  = [5.458472e-02 1.000392e-02 3.480595e-03 1.832472e-03];
L2_p2  = [5.152661e-02 9.918104e-03 2.935698e-03 7.070890e-04];
L2_p3  = [5.123901e-02 1.090484e-02 2.824788e-03 7.134444e-04];
y_vals1= 5*h;
y_vals2= 0.1*h.^2;

figure(1)
clf
hold on
plot(h,L2_pm1,'--*r','LineWidth', 2);
plot(h, L2_p1,'--*g','LineWidth', 2);
plot(h, L2_p2,'--*c','LineWidth', 2);
plot(h, L2_p3,'--*b','LineWidth', 2);
plot(h,y_vals1,'-m','LineWidth', 2);
plot(h,y_vals2,'-k','LineWidth', 2);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('h','Interpreter','latex'); ylabel('$L_2$','Interpreter','latex');
legend('RBF', 'RBF+p=1', 'RBF+p=2', 'RBF+p=3', 'slope=1', 'slope=2','Interpreter','latex');
set(gcf,'color','w');
set(gca, 'FontSize', 18);

%% problem parameters
% near body mesh 2 fringe gap: bottom: 7.362500e-01, right: 6.637500e-01, top: 6.637500e-01, left:7.362500e-01 

% ov_info = containers.Map({ 'num grids', 'mesh1 donor', 'mesh2 donor', 'mandatory frng', 'overlap', 'donor grid', ...
%                            'intrp radius', 'intrp type', 'intrp shape', 'shape param', 'poly order'}, ...
%                          { 2, 2, 1, 2, [3*h2(1), 3*h2(1), 3*h2(1), 3*h2(1)], "radial", ...
%                            2.5*max(h2), "rbf", "gaussian", 1.0, 1 });

% ov_info = containers.Map({ 'num grids', 'mesh1 donor', 'mesh2 donor', 'mandatory frng', 'overlap', 'donor grid', ...
%                            'intrp radius', 'intrp type', 'intrp shape', 'shape param', 'poly order'}, ...
%                          { 2, 2, 1, 3, [7*h2(1), 6*h2(1), 6*h2(1), 7*h2(1)], "radial", ...
%                            2.5*max(h2), "rbf", "gaussian", 1.0, 1 });

% ov_info = containers.Map({ 'num grids', 'mesh1 donor', 'mesh2 donor', 'mandatory frng', 'overlap', 'donor grid', ...
%                            'intrp radius', 'intrp type', 'intrp shape', 'shape param', 'poly order'}, ...
%                          { 2, 2, 1, 5, [14*h2(1), 13*h2(1), 13*h2(1), 14*h2(1)], "radial", ...
%                            2.5*max(h2), "rbf", "gaussian", 1.0, 1 });

% ov_info = containers.Map({ 'num grids', 'mesh1 donor', 'mesh2 donor', 'mandatory frng', 'overlap', 'donor grid', ...
%                            'intrp radius', 'intrp type', 'intrp shape', 'shape param', 'poly order'}, ...
%                          { 2, 2, 1, 9, [29*h2(1), 26*h2(1), 26*h2(1), 29*h2(1)], "radial", ...
%                            2.5*max(h2), "rbf", "gaussian", 1.0, 1 });
                       
del    = 0.2;
h      = [del del/2 del/4 del/8];
L2_pm1 = [5.473572e-02 1.234375e-02 2.959796e-03 7.257568e-04];
L2_p1  = [5.458472e-02 1.240483e-02 2.983622e-03 7.620311e-04];
y_vals1= 5*h;
y_vals2= 0.1*h.^2;

figure(2)
clf
hold on
plot(h, L2_pm1,'--*r','LineWidth', 2);
plot(h, L2_p1,'--*g','LineWidth', 2);
plot(h,y_vals1,'-m','LineWidth', 2);
plot(h,y_vals2,'-k','LineWidth', 2);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('h','Interpreter','latex'); ylabel('$L_2$','Interpreter','latex');
legend('RBF','RBF+p=1','slope=1', 'slope=2','Interpreter','latex');
set(gcf,'color','w');
set(gca, 'FontSize', 18);

%% problem parameters
% near body mesh 2 fringe gap: bottom: 7.362500e-01, right: 6.637500e-01, top: 6.637500e-01, left:7.362500e-01 
%
% ov_info = containers.Map({ 'num grids', 'mesh1 donor', 'mesh2 donor', 'mandatory frng', 'overlap', 'donor grid', ...
%                            'intrp radius', 'intrp type', 'intrp shape', 'shape param', 'poly order'}, ...
%                          { 2, 2, 1, 9, [29*h2(1), 26*h2(1), 26*h2(1), 29*h2(1)], "radial", ...
%                            2.5*max(h2), "rbf", "gaussian", 1.0, 1 });
                       
nnodes = [18,37,64,92,130,176,222,280,340,416,488,565,652,746,854,956,1070,1188,1313,1450,1580,1734,1886,2040,2200];
time   = [0.012,0.014,0.015,0.016,0.019,0.023,0.026,0.033,0.041,0.049,0.06,0.07,0.085,0.1,0.13,0.16,0.2,0.24,0.27,0.33,0.38,0.47,0.54,0.64,0.76];

figure(3)
clf
hold on
plot(nnodes, time/time(1),'--+b','LineWidth', 2);
xlabel('number of nodes','Interpreter','latex'); ylabel('time relative to 1st point','Interpreter','latex');
title('RBF+p=2','Interpreter','latex');
set(gcf,'color','w');
set(gca, 'FontSize', 18);
