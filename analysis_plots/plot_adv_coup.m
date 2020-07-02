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
plot(h,y_vals,  '-k','LineWidth', 2);
plot(h,L2,'--*r','LineWidth', 2,'MarkerSize',20);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('$h$','Interpreter','latex'); ylabel('$L_2$','Interpreter','latex');
legend('slope=2', 'linear Lagrange', 'Interpreter','latex');
set(gcf,'color','w');
set(gca, 'FontSize', 24);

%% periodic bc with non-constant fringe gap

% % boundary condition map - bottom, right, top, left
% bc1 = containers.Map({   'bottom',      'right',       'top',      'left'}, ...
%                      {  "per_mas",  "per_slave", "per_slave",   "per_mas"});

del    = 0.2;
h      = [del del/2 del/4 del/8];
L2_p1  = [1.4619375017366851e-01, 3.9809634740986186e-02, 1.0413186214473841e-02, 2.7245658567234335e-03];
L2_p2  = [1.4298573344364082e-01, 3.8767082590421741e-02, 9.9331847446796715e-03, 2.4999861096006875e-03];
L2_p3  = [1.4435545957188770e-01, 3.8745472213715160e-02, 9.9209346084406288e-03, 2.5014280117732022e-03];
y_vals = 25*h.^2;

figure(3)
clf
hold on
plot(h,y_vals,'-k','LineWidth', 2);
plot(h,L2_p1,'--sr','LineWidth', 2,'MarkerSize',20);
plot(h,L2_p2,'--^g','LineWidth', 2,'MarkerSize',20);
plot(h,L2_p3,'--ob','LineWidth', 2,'MarkerSize',20);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('$h$','Interpreter','latex'); ylabel('$L_2$','Interpreter','latex');
legend('slope=2', 'linear Lagrange', 'quadratic Lagrange', 'cubic Lagrange', 'Interpreter','latex');
set(gcf,'color','w');
set(gca, 'FontSize', 24);

%% periodic bc with non-constant fringe gap using RBF interpolation

% ov_info = containers.Map({ 'num grids', 'mesh1 donor', 'mesh2 donor', 'mandatory frng', 'overlap', 'donor grid', ...
%                            'intrp radius', 'intrp type', 'intrp shape', 'shape param', 'poly order', 'solve type', 'fringe update'}, ...
%                          { 2, 2, 1, 2, [3*h2(1), 3*h2(1), 3*h2(1), 3*h2(1)], "radial", ...
%                            2.5*max(h2), "rbf", "gaussian", 1.0, 1, "coupled", "direct"});

del    = 0.2;
h      = [del del/2 del/4 del/8];
L2_p1  = [1.4462064568919664e-01, 3.8901167384080330e-02, 9.9969732160264703e-03, 2.5302788366218998e-03];
L2_p2  = [1.4405951192382591e-01, 3.8751092903826695e-02, 9.9237695440348075e-03, 2.5010589342802643e-03];
L2_p3  = [1.4440976153230486e-01, 3.8747904430549381e-02, 9.9213383512770890e-03, 2.5014466551637783e-03];
y_vals = 25*h.^2;

figure(4)
clf
hold on
plot(h,y_vals,'-k','LineWidth', 2);
plot(h,L2_p1,'--sr','LineWidth', 2,'MarkerSize',20);
plot(h,L2_p2,'--^g','LineWidth', 2,'MarkerSize',20);
plot(h,L2_p3,'--ob','LineWidth', 2,'MarkerSize',20);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('$h$','Interpreter','latex'); ylabel('$L_2$','Interpreter','latex');
legend('slope=2', 'RBF($C^1$)', 'RBF($C^2$)', 'RBF($C^3$)', 'Interpreter','latex');
set(gcf,'color','w');
set(gca, 'FontSize', 24);

% fig = figure(3);
% set(fig,'Units','Inches');
% pos = get(fig,'Position');
% set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(fig,'temp_name','-dpdf','-r0')
