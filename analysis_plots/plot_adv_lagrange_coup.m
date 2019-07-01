clc; close all; clear all;

% %% problem properties
% pp = containers.Map({                'prblm', 'dof per node', 'velocity'}, ...
%                      { "unsteady scalar adv",              1,      [1,1]});
%                  
% %% background mesh
% box1 = [-2,2; -2,2];
% h1   = [ 0.2,  0.2]/2;
% 
% % boundary condition map - bottom, right, top, left
% bc1 = containers.Map({   'bottom',      'right',       'top',      'left'}, ...
%                      {"dirichlet",  "dirichlet", "dirichlet", "dirichlet"});
% 
% mesh1 = containers.Map({'dim', 'size', 'bc'}, ...
%                        { box1,     h1,  bc1});
%                    
% %% near body mesh
% box2 = [-1,1; -1,1];
% h2   = [ 0.2,  0.2]/2;
% 
% % boundary condition map - bottom, right, top, left
% bc2 = containers.Map({ 'bottom',   'right',     'top',    'left'}, ...
%                      {"overset", "overset", "overset", "overset"});
% 
% mesh2 = containers.Map({'dim', 'size', 'bc'}, ...
%                        { box2,     h2,  bc2});
%                    
% %% implicit hole cutting options -- apply to both grids
% ov_info = containers.Map({ 'num grids', 'background grid', 'mesh1 donor', 'mesh2 donor', 'mandatory frng', 'overlap', 'donor grid', 'intrp order'}, ...
%                          {  2, 1, 2, 1, 2, [4*0.1, 4*0.1, 4*0.1, 4*0.1], "tensor", 1 });
% 
% %% time step and linear solve parameters
% dt = min(h2./pp('velocity'));
% 
% time_sol_info = containers.Map({'init time', 'total time', 'time step', 'BDF order'}, ...
%                              {          0.0,          2.0,          dt,           2} );
% 
% lin_sol_info = containers.Map({'Newton steps', 'residual tolerance'},
%                               {            10,                1e-09});

del    = 0.2;
h      = [del del/2 del/4 del/8];
L2     = [1.3279249867553922e-01, 3.6550624167650321e-02, 9.4535063280176621e-03, 2.3938773619278410e-03];
y_vals = 25*h.^2;

figure(1)
clf
hold on
plot(h,L2,'--*r','LineWidth', 2);
plot(h,y_vals,  '-k','LineWidth', 2);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('h','Interpreter','latex'); ylabel('$L_2$','Interpreter','latex');
title('coincident meshes','Interpreter','latex');
legend('overset p=1', 'slope=2','Interpreter','latex');
set(gcf,'color','w');
set(gca, 'FontSize', 18);

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
L2     = [1.3563407088012361e-01, 3.7807655933261207e-02, 9.9860649996785372e-03, 2.6282128810528730e-03];
y_vals = 25*h.^2;

figure(2)
clf
hold on
plot(h,L2,'--*r','LineWidth', 2);
plot(h,y_vals,  '-k','LineWidth', 2);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('h','Interpreter','latex'); ylabel('$L_2$','Interpreter','latex');
legend('overset p=1', 'slope=2','Interpreter','latex');
set(gcf,'color','w');
set(gca, 'FontSize', 18);
