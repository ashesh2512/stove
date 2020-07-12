clc; close all; clear all;

% %% problem properties
% pp = containers.Map({              'prblm', 'dof per node', 'conductivity', 'velocity', 'mean'}, ...
%                      { "unsteady adv diff",              1,              0,      [1,1],    1.0});
%                  
% %% background mesh
% box1 = [-2,2; -2,2];
% h1   = [ 0.2,  0.2]/2;
% 
% % boundary condition map - bottom, right, top, left
% bc1 = containers.Map({   'bottom',      'right',       'top',      'left'}, ...
%                      {  "per_mas",  "per_slave", "per_slave",   "per_mas"});
%                  
% mesh1 = containers.Map({'dim', 'size', 'bc'}, ...
%                        { box1,     h1,  bc1});
%                    
% %% near body mesh
% box2 = [-1.13625,0.86375; -1.13625,0.86375];
% h2   = [             0.2,              0.2]/2;
% 
% % boundary condition map - bottom, right, top, left
% bc2 = containers.Map({ 'bottom',   'right',     'top',    'left'}, ...
%                      {"overset", "overset", "overset", "overset"});
% 
% mesh2 = containers.Map({'dim', 'size', 'bc'}, ...
%                        { box2,     h2,  bc2});
%
% ov_info = containers.Map({ 'num grids', 'mesh1 donor', 'mesh2 donor', 'mandatory frng', 'overlap', 'donor grid', ...
%                            'intrp radius', 'intrp type', 'intrp shape', 'shape param', 'poly order', 'solve type', 'fringe update'}, ...
%                          { 2, 2, 1, 2, [3*h2(1), 3*h2(1), 3*h2(1), 3*h2(1)], "radial", ...
%                            2.5*max(h2), "rbf", "gaussian", 1.0, 1, "coupled", "direct"});
% 
% %% time step and linear olve parameters
% dt = min(h2./pp('velocity'));
% 
% time_sol_info = containers.Map({'init time', 'total time', 'time step', 'BDF order'}, ...
%                                {        0.0,        1.999,          dt,           2});
% 
% lin_sol_info = containers.Map({'type', 'Newton steps', 'Newton tolerance', 'GMRES iterations', 'GMRES tolerance', ...
%                                'decoupled loops', 'solution tolerance', 'dsol history', 'dsol tolerance'}, ...
%                               {"direct", 3, 1e-09, 1, 1e-6, 120, 1e-1, 10, 1e-15});

del    = 0.2;
h      = [del del/2 del/4 del/8];

L2_C1_pe_inf  = [1.4462064568919650e-01, 3.8901167384079831e-02, 9.9969732160256689e-03 2.5302788366189538e-03,];
L2_C1_pe_1000 = [1.4385469995798469e-01, 3.8683397111696438e-02, 9.9351917965267659e-03, 2.5117182818728281e-03];
L2_C1_pe_100 = [1.3715388499458830e-01, 3.6808030977914657e-02, 9.4376870485671348e-03, 2.3882150872071295e-03];
L2_C1_pe_10 = [8.5233671838704256e-02, 2.2375113931202752e-02, 5.6545916098210327e-03, 1.4189241463203163e-03];

L2_C1_pe_1  = [2.1921410575299486e-03, 3.5310758549192964e-04, 1.6138301015939584e-04, 1.2181816113133137e-04];
L2_C2_pe_1  = [2.1739879491656034e-03, 3.2841466860091879e-04, 6.4658287739179044e-05, 1.4660582742780726e-05];
L2_C3_pe_1  = [2.1655535321836662e-03, 3.2449544864477579e-04, 6.4820292090176325e-05, 1.4621554609798473e-05];

L2_C1_pe_p10  = [2.2472932811554120e-04, 3.4169988176370481e-05, 9.5325415495192410e-05, 8.1318862968400850e-05];
L2_C2_pe_p10  = [1.8860734675442517e-04, 3.6126313521504175e-05, 5.3606163684837370e-07, 6.9658680943685813e-07];
L2_C3_pe_p10  = [4.8345298332457586e-05, 1.8203051870485835e-06, 2.4460283266438742e-06, 4.2645817898455545e-07];

h_more      = [del del/2 del/4 del/8 del/16];
L2_C1_pe_0  = [2.0132630043953880e-04, 2.9092760191521902e-05, 8.8983158302284894e-05, 7.5976056795078629e-05];
L2_C2_pe_0  = [1.8166748083944174e-04, 3.5651015544360451e-05, 3.5658491486951507e-07, 6.6251707792438454e-07, 2.4375326290910114e-08];

y_vals1 = 100*h.^2;
y_vals2 = 1e-5*h.^2;


figure(2)
clf
hold on

plot(h,y_vals1,  '-k','LineWidth', 2);
plot(h,L2_C1_pe_1,'--*r','LineWidth', 2,'MarkerSize',20);
plot(h,L2_C2_pe_1,'--^g','LineWidth', 2,'MarkerSize',20);
plot(h,L2_C3_pe_1,'--ob','LineWidth', 2,'MarkerSize',20);

plot(h,y_vals2,  '-k','LineWidth', 2);
plot(h,L2_C1_pe_100,'--*r','LineWidth', 2,'MarkerSize',20);
plot(h,L2_C1_pe_10,'--*r','LineWidth', 2,'MarkerSize',20);
plot(h,L2_C1_pe_p10,'--*r','LineWidth', 2,'MarkerSize',20);
plot(h,L2_C2_pe_p10,'--^g','LineWidth', 2,'MarkerSize',20);
plot(h,L2_C3_pe_p10,'--ob','LineWidth', 2,'MarkerSize',20);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('$h$','Interpreter','latex'); ylabel('$L_2$','Interpreter','latex');
legend('slope=2', 'RBF($C^1$)', 'RBF($C^2$)', 'RBF($C^3$)', 'Interpreter','latex');
set(gcf,'color','w');
set(gca, 'FontSize', 24);

th = text(0.1,1e-1,0,'$Pe/L = 100$');
set(th,'FontSize',20,'BackgroundColor','none','EdgeColor','none','Interpreter','latex');
th = text(0.1,1e-2,0,'$Pe/L = 10$');
set(th,'FontSize',20,'BackgroundColor','none','EdgeColor','none','Interpreter','latex');
th = text(0.1,1e-3,0,'$Pe/L = 1$');
set(th,'FontSize',20,'BackgroundColor','none','EdgeColor','none','Interpreter','latex');
th = text(0.1,2e-5,0,'$Pe/L = 0.1$');
set(th,'FontSize',20,'BackgroundColor','none','EdgeColor','none','Interpreter','latex');


