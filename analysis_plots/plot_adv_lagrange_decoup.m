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
% box2 = [-1.13625,0.86375; -1.13625,0.86375];
% h2   = [             0.2,              0.2];
% 
% % boundary condition map - bottom, right, top, left
% bc2 = containers.Map({ 'bottom',   'right',     'top',    'left'}, ...
%                      {"overset", "overset", "overset", "overset"});
% 
% mesh2 = containers.Map({'dim', 'size', 'bc'}, ...
%                        { box2,     h2,  bc2});
%                    
% %% implicit hole cutting options -- apply to both grids
% ov_info = containers.Map({ 'num grids', 'background grid', 'mesh1 donor', 'mesh2 donor', 'mandatory frng', 'overlap', 'donor grid', 'intrp order', 'solve type'}, ...
%                          {  2, 1, 2, 1, 2, [3*h2(1), 3*h2(1), 3*h2(1), 3*h2(1)], "tensor", 1, "decoupled" });
% 
% %% time step and linear solve parameters
% dt = min(h2./pp('velocity'));
% 
% time_sol_info = containers.Map({'init time', 'total time', 'time step', 'BDF order'}, ...
%                              {          0.0,          2.0,          dt,           2} );
% 
% lin_sol_info = containers.Map({'Newton steps', 'residual tolerance', 'decoupled loops', 'solution tolerance'}, ...
%                               {             3,                1e-06,               120,                 1e-1});

del     = 0.2;
h       = [del del/2 del/4 del/8];
L2      = [1.3563407088012361e-01, 3.7807655933261207e-02, 9.9860649996785372e-03, 2.6282128810528730e-03];
L2_1em1 = [1.3563575177060619e-01, 3.7829257539133951e-02, 9.9844824580420759e-03, 2.5954002338079477e-03];
L2_1em4 = [1.3563402647239917e-01, 3.7807621123918798e-02, 9.9859756176621326e-03, 2.6280764437477504e-03];
y_vals  = 25*h.^2;

figure(1)
clf
hold on
plot(h,y_vals,  '-k','LineWidth', 2);
plot(h,L2,'--*r','LineWidth', 2,'MarkerSize',20);
plot(h,L2_1em1,'--sb','LineWidth', 2,'MarkerSize',20);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('$h$','Interpreter','latex'); ylabel('$L_2$','Interpreter','latex');
legend('slope=2','tightly coupled', '2 loosely coupled iterations','Interpreter','latex');
set(gcf,'color','w');
set(gca, 'FontSize', 24);

%% periodic bc with non-constant fringe gap

% % boundary condition map - bottom, right, top, left
% bc1 = containers.Map({   'bottom',      'right',       'top',      'left'}, ...
%                      {  "per_mas",  "per_slave", "per_slave",   "per_mas"});

del     = 0.2;
h       = [del del/2 del/4 del/8];
L2      = [1.4619375017366851e-01, 3.9809634740986186e-02, 1.0413186214473841e-02, 2.7245658567234335e-03];
L2_1em1 = [1.4619528859468175e-01, 3.9830139066427690e-02, 1.0411670122109485e-02, 2.6929274388159362e-03];
y_vals  = 25*h.^2;

figure(2)
clf
hold on
plot(h,y_vals,  '-k','LineWidth', 2);
plot(h,L2,'--*r','LineWidth', 2,'MarkerSize',20);
plot(h,L2_1em1,'--sb','LineWidth', 2,'MarkerSize',20);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('$h$','Interpreter','latex'); ylabel('$L_2$','Interpreter','latex');
legend('slope=2','OSS', 'OAS, $2OC$','Interpreter','latex');
set(gcf,'color','w');
set(gca, 'FontSize', 20);

% fig = figure(2);
% set(fig,'Units','Inches');
% pos = get(fig,'Position');
% set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(fig,'temp_name','-dpdf','-r0')
