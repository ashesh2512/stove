clc; clear all; close all;

addpath(genpath('../src')); % add source folder

tol = 1e-14; %set regression test tolerance

total_pass = 0; % counter for tests passed
total_fail = 0; % counter for tests failed

%% test heat equation on coupled meshes using lagrange interpolation
[num_pass,num_fail] = heat_coup_lag(tol);

total_pass = total_pass+num_pass;
total_fail = total_fail+num_fail;

%% test heat equation on decoupled meshes using lagrange interpolation
[num_pass,num_fail] = heat_decoup_lag(tol);

total_pass = total_pass+num_pass;
total_fail = total_fail+num_fail;

%% test heat equation on coupled meshes using rbf interpolation
[num_pass,num_fail] = heat_coup_rbf(tol);

total_pass = total_pass+num_pass;
total_fail = total_fail+num_fail;

%% test advection equation on coupled meshes using lagrange interpolation
[num_pass,num_fail] = adv_coup_lag(tol);

total_pass = total_pass+num_pass;
total_fail = total_fail+num_fail;

%% test advection equation on coupled meshes using rbf interpolation
[num_pass,num_fail] = adv_coup_rbf(tol);

total_pass = total_pass+num_pass;
total_fail = total_fail+num_fail;

%% test advection equation on decoupled meshes using lagrange interpolation
[num_pass,num_fail] = adv_decoup_lag(tol);

total_pass = total_pass+num_pass;
total_fail = total_fail+num_fail;

%% test advection diffusion equation on decoupled meshes using lagrange interpolation
[num_pass,num_fail] = adv_diff_coup(tol);

total_pass = total_pass+num_pass;
total_fail = total_fail+num_fail;

%% tally status of all tests

fprintf('\n');
fprintf('Number of tests passed = %u', total_pass);
fprintf('\n');
fprintf('Number of tests failed = %u', total_fail);
fprintf('\n');