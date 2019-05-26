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
% 
% % implicit hole cutting options -- apply to both grids
% % mandatory firnge - number of points inward from the overset boundary to
% %                    be marked as fringe
% % overlap          - minimum distance between innermost fringe points on 
% %                    each grid this value can vary from boundary to 
% %                    boundary
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
plot(h,L2,'--*r','LineWidth', 2);
plot(h,y_vals,  '-k','LineWidth', 2);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('h','FontSize',18); ylabel('L_2','FontSize',18);
title('coincident meshes','FontSize',18);
legend('overset p=1', 'slope=2');
set(gcf,'color','w');
set(gca, 'FontSize', 18);


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
% near body mesh
% box2 = [-1.13625,0.86375; -1.13625,0.86375];
% h2   = [             0.2,              0.2];
% 
% % boundary condition map - bottom, right, top, left
% bc2 = containers.Map({ 'bottom',   'right',     'top',    'left'}, ...
%                      {"overset", "overset", "overset", "overset"});
% 
% % implicit hole cutting options -- apply to both grids
% % mandatory firnge - number of points inward from the overset boundary to
% %                    be marked as fringe
% % overlap          - minimum distance between innermost fringe points on 
% %                    each grid this value can vary from boundary to 
% %                    boundary
% ov_info = containers.Map({ 'num grids', 'mesh1 donor', 'mesh2 donor', 'mandatory frng', 'overlap', 'donor_grid', 'intrp_type', 'intrp_order'}, ...
%                          {  2, 2, 1, 2, [4*0.1, 4*0.1, 4*0.1, 4*0.1], "tensor", "lagrange", 1 });
% 
% % convergence parameters
% N_iters    =    10; % Maximum number of Newton steps
% resnrmdrop = 1e-09; % Newton convergence criteria

del    = 0.2;
h      = [del del/2 del/4 del/8];
L2     = [4.598126e-02 1.113857e-02 2.782573e-03 6.979062e-04];
L2_p1  = [6.526605e-02 9.348554e-03 6.497774e-03 5.626061e-03];
L2_p2  = [4.684465e-02 8.873482e-03 3.253635e-03 6.685104e-04];
L2_p3  = [4.962786e-02 1.088594e-02 2.814165e-03 7.110412e-04];
y_vals1= 5*h;
y_vals2= 0.1*h.^2;

figure(2)
clf
hold on
plot(h,L2_p1,'--*r','LineWidth', 2);
plot(h,L2_p2,'--*g','LineWidth', 2);
plot(h,L2_p3,'--*b','LineWidth', 2);
plot(h,y_vals1,'-m','LineWidth', 2);
plot(h,y_vals2,'-k','LineWidth', 2);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('h','FontSize',18); ylabel('L_2','FontSize',18);
legend('overset p=1', 'overset p=2', 'overset p=3', 'slope=1', 'slope=2');
set(gcf,'color','w');
set(gca, 'FontSize', 18);