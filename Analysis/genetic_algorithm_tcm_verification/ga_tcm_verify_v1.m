clear;
clear global;
close all; 
clc;
format longg;
addpath(genpath('../../'));
load('nri_inclined_x0_xopt.mat');

%% x0

[stm_i0, stt_i0, x_i_f0, x_t0, stm_t0, stt_t_i, t0, t_s0, stm_t_i0]  = createStateStmSttHistory(simparams.x0, simparams);
[tcm_time0, tcm_idx0, min_tcm_dv0, ~, ~, tcm_dv_each0] = opt_multiple_tcm(simparams.x0, t0, t_s0, stm_t0, simparams); % inputs: x, t, t_s, stm_t, stm_t_i, simparams


%% ga
x0 = reshape(simparams.x0,simparams.m,simparams.n);

% Find the indices between the first two maneuvers
t_dv1 = 0;
t_dv2 = sum(x0(7,1:simparams.maneuverSegments(2)-1));
t_dv3 = sum(x0(7,1:simparams.maneuverSegments(3)-1));

idx_dv1 = 1;
idx_dv2 = find(t0 == t_dv2);
idx_dv3 = find(t0 == t_dv3);

% # of TCMs in traj portion 1
tcms_in_p1 = tcm_idx0 < idx_dv2;
num_tcms_p1 = sum(tcms_in_p1);

% # of TCMs in traj portion 2
tcms_in_p2 = tcm_idx0 >= idx_dv2;
num_tcms_p2 = sum(tcms_in_p2);

int_con1 = (1:num_tcms_p1)';

% Bounds
lb = idx_dv1:num_tcms_p1;
ub = idx_dv2-num_tcms_p1+1:idx_dv2;

% Linear inequality constraints
A = eye(num_tcms_p1 - 1, num_tcms_p1);
Im = eye(num_tcms_p1 - 1);
AmI = zeros(num_tcms_p1 - 1, num_tcms_p1);
AmI(:, 2:end) = Im;

A = A-AmI;

b = -1 * ones(num_tcms_p1 - 1, 1);


% NonLinear Inequality constraints
% (sequential) 
% [cin, ceq] = sequential_constraint(tcm_idx0(1:4))

% Objective function
% testTCMdv = calc_tcm_dv(tcm_idx0, x0, t0, t_s0, stm_t0, simparams.P_initial, idx_dv1, idx_dv2, simparams)

%% GA Call and options

% options = optimoptions(@ga,...
%             		'ConstraintTolerance',1e-10,...
%             		'FunctionTolerance',1e-10,...
%             		'MaxTime',36000,...
%             		'PopulationSize',100,...
%            		    'MaxGenerations',500 ,...
%            		    'Display','iter',...
%             		'UseParallel',true,...
%             		'PlotFcn',{@gaplotstopping,...
%                        		@gaplotscores,...
%                        		@gaplotscorediversity,...
%                        		@gaplotbestf,...
%                        		@gaplotbestindiv,...
%                       		@gaplotexpectation});

options = optimoptions(@ga,...
            		'ConstraintTolerance',1e-12,...
            		'FunctionTolerance',1e-12,...
            		'MaxTime',36000,...
            		'PopulationSize',100,...
           		    'MaxGenerations',500 ,...
           		    'Display','iter',...
            		'UseParallel',true);




[tcm_idx_ga, fval] = ga(@(tcm_idx_test)calc_tcm_dv(tcm_idx_test, x0, t0, t_s0, stm_t0, simparams.P_initial, idx_dv1, idx_dv2, simparams), num_tcms_p1, A, b, [], [], lb, ub, @(tcm_idx_test)sequential_constraint(tcm_idx_test, x0, t0, t_s0, stm_t0, simparams.P_initial, idx_dv1, idx_dv2, simparams), int_con1, options)












% tcm_idx_ga = [7, 474, 1856, 2287];

sequential_constraint(tcm_idx_ga, x0, t0, t_s0, stm_t0, simparams.P_initial, idx_dv1, idx_dv2, simparams)


%% x_opt

% [stm_i, stt_i, x_i_f, x_t, stm_t, stt_t_i, t, t_s, stm_t_i]  = createStateStmSttHistory(x_opt, simparams);
% [tcm_time,tcm_idx,min_tcm_dv] = opt_multiple_tcm(x_opt, t, t_s, stm_t, simparams);








%% constraint function



%% objective function - calc_tcm_dv


















