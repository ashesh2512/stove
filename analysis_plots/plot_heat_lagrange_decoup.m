clc; close all; clear all;

%% problem properties
% pp = containers.Map({            'prblm', 'dof per node', 'conductivity'}, ...
%                      { "steady heat MMS",              1,              1});
%                  
% % background mesh
% box1 = [-2,2; -2,2];
% h1   = [ 0.2,  0.2];
% 
% % boundary condition map - bottom, right, top, left
% bc1 = containers.Map({   'bottom',     'right',       'top',      'left'}, ...
%                      {"dirichlet", "dirichlet", "dirichlet", "dirichlet"});
% 
% mesh1 = containers.Map({'dim', 'size', 'bc'}, ...
%                        { box1,     h1,  bc1});
%                    
% % near body mesh
% box2 = [-1,1; -1,1];
% h2   = [ 0.2,  0.2];
% 
% % boundary condition map - bottom, right, top, left
% bc2 = containers.Map({ 'bottom',   'right',     'top',    'left'}, ...
%                      {"overset", "overset", "overset", "overset"});
% 
% mesh2 = containers.Map({'dim', 'size', 'bc'}, ...
%                        { box2,     h2,  bc2});
%                    
% % implicit hole cutting options -- apply to both grids
% ov_info = containers.Map({ 'num grids', 'mesh1 donor', 'mesh2 donor', 'mandatory frng', 'overlap', 'donor grid', 'intrp order', 'solve type'}, ...
%                          {  2, 2, 1, 2, [4*0.1, 4*0.1, 4*0.1, 4*0.1], "tensor", 1, "decoupled" });
% 
% time step and linear solve parameters
% 
% time_sol_info = containers.Map({'init time', 'total time', 'time step', 'BDF order'}, ...
%                              {          0.0,          2.0,         100,           2});
% 
% lin_sol_info = containers.Map({'Newton steps', 'residual tolerance', 'decoupled loops', 'solution tolerance'}, ...
%                               {            10,                1e-09,                10,                1e-03});

del = 0.2;
h   = [del del/2 del/4 del/8];
L2     = [4.598126e-02 1.113857e-02 2.782573e-03 6.979062e-04];
L2_tol_1em1 = [5.3772858819218880e-02, 3.7635260003548165e-02, 4.0748274281320176e-02, 4.1376668209188297e-02];
L2_tol_1em2 = [4.4219744584391485e-02, 1.1391949504985075e-02, 5.0662117370762491e-03, 4.4193950330099728e-03];
L2_tol_1em3 = [4.5862781242750561e-02, 1.1114255449853644e-02, 2.7992868742053936e-03, 7.8918339502312234e-04];
y_vals = 15*h.^2;

figure(1)
clf
hold on
plot(h,L2,'--*m','LineWidth', 2);
plot(h,L2_tol_1em1,'--*r','LineWidth', 2);
plot(h,L2_tol_1em2,'--*b','LineWidth', 2);
plot(h,L2_tol_1em3,'--*c','LineWidth', 2);
plot(h,y_vals,  '-k','LineWidth', 2);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('h','Interpreter','latex'); ylabel('$L_2$','Interpreter','latex');
title('decoupled solve on coincident meshes','Interpreter','latex');
legend('coupled','tol = 1e-01', 'tol = 1e-02', 'tol = 1e-03', 'slope=2','Interpreter','latex');
set(gcf,'color','w');
set(gca, 'FontSize', 18);

L2_hp2   = [0.0826761802911245,0.0639114897012161,0.0537728588192189,0.0501767610782145,0.0435758693746887,0.0451531260523086,0.0415329414415355,0.0441400252647475,0.0422018749942269,0.0443680070834562,0.0433060674015262,0.0448200177698567,0.0442197445843915,0.0452086993875077,0.0448605883484789,0.0454866184185041,0.0452810584553684,0.0456709927414977,0.0455481509387044,0.0457888835705408,0.0457149119248654,0.0458627812427506,0.0458180244553727,0.0459085836316680,0.0458814231130731];
L2_hp1   = [0.132030015668975,0.0990344473121247,0.0769628088008273,0.0606744323069027,0.0471610455296539,0.0376352600035482,0.0295577904270041,0.0242365337414505,0.0196127818785449,0.0169246524401160,0.0145007041394406,0.0133649201321402,0.0122367800290822,0.0118695275435007,0.0113919495049851,0.0113243219006112,0.0111298448316186,0.0111518463286386,0.0110721373070086,0.0111093573180794,0.0110755471109009,0.0111069884837612,0.0110918946537978,0.0111142554498536,0.0111071017476446];
L2_hp05  = [0.148370809508374,0.110563204602723,0.0856806229269238,0.0670967108630976,0.0520091730223828,0.0407482742813202,0.0316211507908779,0.0247957652774845,0.0192989411872312,0.0151783775511657,0.0119016293343352,0.00944309559750964,0.00753741766953121,0.00611477342032622,0.00506621173707625,0.00429592954720232,0.00377550207709237,0.00340041505682849,0.00317597263508641,0.00301227487797820,0.00292801109160964,0.00286128780704800,0.00283350137360959,0.00280685691415313,0.00279928687420539];
L2_hp025 = [0.153687390386786,0.112788744791183,0.0880095223436267,0.0682891571316498,0.0533136545199633,0.0413766682091883,0.0323101390016447,0.0250727368561148,0.0195884597965757,0.0151988729295068,0.0118860758621716,0.00922332540537432,0.00722792174741009,0.00561381083707289,0.00441939503300997,0.00344449610548616,0.00273962173832429,0.00215783949917802,0.00175508786576935,0.00141886059973672,0.00120386654697656,0.00102208767904537,0.000920070745247008,0.000830771051943650,0.000789183395023122];

iterations = 1:length(L2_hp025);
oneD_arr = ones(length(L2_hp025));
y_vals = h.^2;

figure(2)
clf
h1 = plot3(iterations,0.2*oneD_arr,L2_hp2,'--*r','LineWidth', 2);
hold on
h2 = plot3(iterations,0.1*oneD_arr,L2_hp1,'--*m','LineWidth', 2);
h3 = plot3(iterations,0.05*oneD_arr,L2_hp05,'--*b','LineWidth', 2);
h4 = plot3(iterations,0.025*oneD_arr,L2_hp025,'--*c','LineWidth', 2);
h5 = plot3([ 5  5  5  5],h,1.1*h.^2,'-k','LineWidth', 2);
plot3([15 15 15 15],h,1.1*h.^2,'-k','LineWidth', 2);
plot3([25 25 25 25],h,1.2*h.^2,'-k','LineWidth', 2);
grid on;
set(gca, 'YScale', 'log', 'ZScale', 'log');
xlabel('decoupled iterations','Interpreter','latex'); ylabel('h','Interpreter','latex'); zlabel('$L_2$','Interpreter','latex');
legend([h1(1),h2(1),h3(1),h4(1),h5(1)], 'h = 0.2','h = 0.1','h = 0.05','h = 0.025','slope=2','Interpreter','latex');
set(gcf,'color','w');
set(gca, 'FontSize', 18);
view(164.8603,23.6139);

% fig = figure(2);
% set(fig,'Units','Inches');
% pos = get(fig,'Position');
% set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(fig,'iterations_error_tradeoff','-dpdf','-r0')

%% problem properties

% box2 = [-1.13625,0.86375; -1.13625,0.86375];

% ov_info = containers.Map({ 'num grids', 'mesh1 donor', 'mesh2 donor', 'mandatory frng', 'overlap', 'donor grid', 'intrp order', 'solve type'}, ...
%                          {  2, 2, 1, 2, [3*h2(1), 3*h2(1), 3*h2(1), 3*h2(1)], "tensor", 2, "decoupled" });


del = 0.2;
h   = [del del/2 del/4 del/8];
L2_p1  = [6.526605e-02 9.348554e-03 6.497774e-03 5.626061e-03];
L2_p2  = [4.684465e-02 8.873482e-03 3.253635e-03 6.685104e-04];
L2_tol_1em3_p1 = [6.5337221723682393e-02, 9.3936384501601099e-03, 6.5838580129703944e-03, 7.9999507815616232e-03];
L2_tol_1em4_p1 = [6.5284671745566167e-02, 9.3419045678414282e-03, 6.5047293687424896e-03, 5.8161874023218412e-03];
L2_tol_1em3_p2 = [4.7027290105483073e-02, 8.8711570598291236e-03, 3.5933359389325823e-03, 3.0210264350805185e-03];
L2_tol_1em4_p2 = [4.6865713243357733e-02, 8.8748808342153582e-03, 3.2637912890110957e-03, 8.2209112619829698e-04];
y_vals = 15*h.^2;


figure(3)
clf
subplot(1,2,1)
hold on
plot(h,L2_p1,'--*m','LineWidth', 2);
plot(h,L2_tol_1em3_p1,'--*r','LineWidth', 2);
plot(h,L2_tol_1em4_p1,'--*b','LineWidth', 2);
plot(h,y_vals,  '-k','LineWidth', 2);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('h','Interpreter','latex'); ylabel('$L_2$','Interpreter','latex');
legend('coupled', 'tol = 1e-03, decoupled', 'tol = 1e-04, decoupled', 'slope=2','Interpreter','latex');
set(gcf,'color','w');
set(gca, 'FontSize', 18);

subplot(1,2,2)
hold on
plot(h,L2_p2,'--*m','LineWidth', 2);
plot(h,L2_tol_1em3_p2,'--*r','LineWidth', 2);
plot(h,L2_tol_1em4_p2,'--*b','LineWidth', 2);
plot(h,y_vals,  '-k','LineWidth', 2);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('h','Interpreter','latex'); ylabel('$L_2$','Interpreter','latex');
legend('coupled', 'tol = 1e-03, decoupled', 'tol = 1e-04, decoupled', 'slope=2','Interpreter','latex');
set(gcf,'color','w');
set(gca, 'FontSize', 18);
