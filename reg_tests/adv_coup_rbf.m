function [test_passed,test_failed] = adv_coup_rbf(tol)

% set counter for passed and failed tests
test_passed = 0;
test_failed = 0;

%% problem properties
pp = containers.Map({                'prblm', 'dof per node', 'velocity'}, ...
                     { "unsteady scalar adv",              1,      [1,1]});
                 
%% background mesh
box1 = [-2,2; -2,2];
h1   = [ 0.2,  0.2]/2;

% boundary condition map - bottom, right, top, left
bc1 = containers.Map({   'bottom',      'right',       'top',      'left'}, ...
                     {  "per_mas",  "per_slave", "per_slave",   "per_mas"});
                 
mesh1 = containers.Map({'dim', 'size', 'bc'}, ...
                       { box1,     h1,  bc1});
                   
%% near body mesh
box2 = [-1.13625,0.86375; -1.13625,0.86375];
h2   = [             0.2,              0.2]/2;

% boundary condition map - bottom, right, top, left
bc2 = containers.Map({ 'bottom',   'right',     'top',    'left'}, ...
                     {"overset", "overset", "overset", "overset"});

mesh2 = containers.Map({'dim', 'size', 'bc'}, ...
                       { box2,     h2,  bc2});
                   
%% implicit hole cutting options -- apply to both grids
% num grids        - number of grids is currently restricted to 2
% mesh1 donor      - Donor mesh for mesh 1; currently this can only be 2
% mesh2 donor      - Donor mesh for mesh 2; currently this can only be 1
% mandatory firnge - number of points inward from the overset boundary to
%                    be marked as fringe
% overlap          - minimum distance between innermost fringe points on 
%                    each grid this value can vary from boundary to 
%                    boundary
% donor grid       - Donor grid type; options are tensor or radial
% intrp radius     - Size of donor grid if radial selected
% intrp type       - Type of interpolation desired. Tensor grid is hard
%                    coded to performe a Lagrange interpolation. Radial 
%                    grid has only RBF as an option right now.
% intrp shape      - For RBF this determined the type of RBF. Gaussian is
%                    only option currently.
% shape param      - Shape parameter for RBF interpolation.
% intrp order      - Order of consistency desired in interpolating
%                    functions. -1 for RBF uses a classical RBF with 0th 
%                    order consistency.
ov_info = containers.Map({ 'num grids', 'mesh1 donor', 'mesh2 donor', 'mandatory frng', 'overlap', 'donor grid', 'intrp order', 'solve type','fringe update'}, ...
                         {  2, 2, 1, 2, [3*h2(1), 3*h2(1), 3*h2(1), 3*h2(1)], "tensor", 3, "coupled", "direct" });

%% time step and linear solve parameters
dt = min(h2./pp('velocity'));

time_sol_info = containers.Map({'init time', 'total time', 'time step', 'BDF order'}, ...
                             {          0.0,         1.99,          dt,           2} );

lin_sol_info = containers.Map({  'type', 'Newton steps', 'Newton tolerance'}, ...
                              {"direct",             10,              1e-09} );

%% debug/display flags
debug_flags = containers.Map({'plot mesh', 'plot hole cut', 'print fringe gap', 'plot sol'}, ...
                             {      false,           false,              false,      false});

%% coupled meshes, lagrange interpolation, p3

fprintf('\n');
fprintf(['Starting test for scalar linear advection equation on coupled meshes ', ...
         'using lagrange interpolation of order 1 ']);
fprintf('\n');
                                        
%create container map for easy addition and removal of variables without
%having to change other input files
inp_container = containers.Map({'problem definition', 'mesh 1', 'mesh 2', 'overset prop', 'time solver prop', 'lin solver prop', 'debug flags'}, ...
                                {pp, mesh1, mesh2, ov_info, time_sol_info, lin_sol_info, debug_flags} );

L2_err = driver(inp_container);
gold   = 3.8745472213715160e-02;

% set pass/fail status
if abs(L2_err - gold) < tol
    test_status = 'passed'; test_passed = test_passed+1;
else
    test_status = 'failed'; test_failed = test_failed+1;
end
fprintf('\n');
fprintf(['Test for advection equation on coupled coincident meshes ', ...
         'using lagrange interpolation of order 1 ', test_status]);
fprintf('\n');
fprintf('\n');

end