clear;
clear global;
close all; 
clc;
format longg;
addpath(genpath('../../'));



% init_fn = 'init_simpar_nri_3dv_TCMsNodes_flybyNotCorrected';
% run(init_fn);
% x = reshape(simparams.x0,simparams.m,simparams.n);



% load('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\20240229_0837.29_3dv_nri_meddx0_flybynotcorrected_det\workspace.mat')
% x=reshape(x_opt, simparams.m, simparams.n);


savename = '';

%% x0

% [stm_i0, stt_i0, x_i_f0, x_t0, stm_t0, stt_t_i, t0, t_s0, stm_t_i0]  = createStateStmSttHistory(simparams.x0, simparams);
% [tcm_time0, tcm_idx0, min_tcm_dv0, ~, ~, tcm_dv_each0] = opt_multiple_tcm(simparams.x0, t0, t_s0, stm_t0, simparams); % inputs: x, t, t_s, stm_t, stm_t_i, simparams



% commented this out when loading a workspace directly below
% traj  = createStateStmSttQdQHistory(x, simparams);
% [deltaV, deltaVs_nom] = calcDeltaV(x, traj.x_i_f, traj.stm_i, simparams);
% 
% [tcm_time, tcm_idx, min_tcm_dv, ~, ~, tcm_dv_each] = opt_multiple_tcm_wQ(x, traj, deltaVs_nom, simparams); % inputs: x, t, t_s, stm_t, stm_t_i, simparams

load('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\20240228_1754.43_3dv_nri_meddx0_flybynotcorrected_robust\workspace.mat')


% load('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\20240229_0837.29_3dv_nri_meddx0_flybynotcorrected_det\workspace.mat')

% load('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\20240229_0837.29_3dv_nri_meddx0_flybynotcorrected_det\workspace.mat')

% load('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\20240222_1727.03_2dv_leo_llo_det\workspace.mat')
x = reshape(x_opt, simparams.m, simparams.n);


[event_times, event_indicator] = define_events_v2(x(:), traj.t, tcm_time, simparams);

event_idx_logical = logical(sum(traj.t'==event_times', 1));    
event_idxs = find(event_idx_logical);
%% ga

% Find the indices between the first two maneuvers
t_dv1 = 0;
t_dv1 = sum(x(7,1:simparams.maneuverSegments(1)-1));
t_dv2 = sum(x(7,1:simparams.maneuverSegments(2)-1));

idx_dv1 = find(traj.t == t_dv1);
idx_dv2 = find(traj.t == t_dv2);


% # of TCMs in traj portion 1
tcms_in_p1 = tcm_idx < idx_dv2;
num_tcms_p1 = sum(tcms_in_p1);

if length(simparams.maneuverSegments)>2
    t_dv3 = sum(x(7,1:simparams.maneuverSegments(3)-1));

    idx_dv3 = find(traj.t == t_dv3);


    % # of TCMs in traj portion 2
    tcms_in_p2 = tcm_idx >= idx_dv2;
    num_tcms_p2 = sum(tcms_in_p2);


end





int_con1 = (1:num_tcms_p1)';

% Bounds
lb = idx_dv1+1:idx_dv1 + num_tcms_p1;
ub = idx_dv2-num_tcms_p1:idx_dv2 - 1;

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



% start_node = 2;
% target_node = simparams.maneuverSegments(2);
% [tcm_idx_ga, fval_ga] = ga(@(tcm_idx_test)calc_tcm_dv_wQ(tcm_idx_test, x, traj, simparams.P_initial, idx_dv1, idx_dv2, start_node, target_node, 0, deltaVs_nom, simparams), num_tcms_p1, A, b, [], [], lb, ub, @(tcm_idx_test)sequential_constraint_wQ(tcm_idx_test, x, traj, simparams.P_initial, idx_dv1, idx_dv2, start_node, target_node, 0, deltaVs_nom, simparams), int_con1, options)
if length(simparams.maneuverSegments) == 2
    tic
    [tcm_idx_ga, fval_ga] = ga(@(tcm_idx_test)calc_tcm_dv_wQ_entireTraj(tcm_idx_test, x, traj, deltaVs_nom, simparams), num_tcms_p1, A, b, [], [], lb, ub, @(tcm_idx_test)sequential_constraint_wQ_entireTraj(tcm_idx_test, x, traj, deltaVs_nom, simparams), int_con1, options)
    toc

    fval_ga*ndVel2kms*3000
    min_tcm_dv*ndVel2kms*3000
end
% Compare cost of GA to opt:
% tcm_idx_seg1 = tcm_idx(tcms_in_p1);
% tcm_dv_seg1 = calc_tcm_dv_wQ(tcm_idx_seg1, x, traj, simparams.P_initial, idx_dv1, idx_dv2, start_node, target_node, 0, deltaVs_nom, simparams)





% tcm_idx_ga = [7, 474, 1856, 2287];

% sequential_constraint(tcm_idx_ga, x0, t0, t_s0, stm_t0, simparams.P_initial, idx_dv1, idx_dv2, simparams)


%% Now try on entire trajectory
% Linear inequality constraints
% A = eye(num_tcms_p1 - 1, num_tcms_p1);
% Im = eye(num_tcms_p1 - 1);
% AmI = zeros(num_tcms_p1 - 1, num_tcms_p1);
% AmI(:, 2:end) = Im;
% 
% A = A-AmI;
% 
% b = -1 * ones(num_tcms_p1 - 1, 1);




if length(simparams.maneuverSegments)>2
    
    % Linear inequality constraints
    A2 = eye(num_tcms_p2 - 1, num_tcms_p2);
    Im2 = eye(num_tcms_p2 - 1);
    AmI2 = zeros(num_tcms_p2 - 1, num_tcms_p2);
    AmI2(:, 2:end) = Im2;
    
    A2 = A2-AmI2;
    
    b2 = -1 * ones(num_tcms_p2 - 1, 1);
    
    
    
    % Bounds
    % lb = idx_dv1+1:idx_dv1 + num_tcms_p1;
    % ub = idx_dv2-num_tcms_p1:idx_dv2 - 1;
    lb2 = idx_dv2 + 1:idx_dv2 + num_tcms_p2;
    ub2 = idx_dv3 - num_tcms_p2:idx_dv3 - 1;
    
    lbc = [lb, lb2];
    ubc = [ub, ub2];
    
    
    % Combined 1 and 2
    num_tcms = length(tcm_time);
    Ac = zeros(num_tcms - 2, num_tcms);
    Ac(1:num_tcms_p1-1, 1:num_tcms_p1) = A;
    Ac(num_tcms_p1:end, num_tcms_p1+1:end) = A2;
    
    
    bc = [b; b2];
    
    
    int_concc = (1:num_tcms)';
    tic
    [tcm_idx_ga_entireTraj, fval_ga_entireTraj] = ga(@(tcm_idx_test)calc_tcm_dv_wQ_entireTraj(tcm_idx_test, x, traj, deltaVs_nom, simparams), num_tcms_p1 + num_tcms_p2, Ac, bc, [], [], lbc, ubc, @(tcm_idx_test)sequential_constraint_wQ_entireTraj(tcm_idx_test, x, traj, deltaVs_nom, simparams), int_concc, options)
    toc

    fval_ga_entireTraj * ndVel2kms*3000
    min_tcm_dv * ndVel2kms*3000
    
    tcm_time_ga = traj.t(tcm_idx_ga_entireTraj)';
    
    Q_k_km1 = calc_Q_events(traj, x, tcm_time_ga, simparams);
    
    [Pga, min_tcm_dv_totalga, tcm_dv_eachga, P_i_minusga, P_i_plusga] = calc_covariance_wQ_tcmdv_v3(x(:), traj, tcm_time_ga, 1, deltaVs_nom, simparams.P_initial, Q_k_km1, simparams);
end

%% constraint function



%% objective function - calc_tcm_dv


%% save output


if 1
    formatOut = 'yyyymmdd_HHMM.SS';
    dateString = datestr(now,formatOut);
    outputPath = strcat('./sims/',dateString,'_',savename);
    mkdir(outputPath);
    save(strcat('./',outputPath,'/workspace.mat'));

    % save main file backup
    FileNameAndLocation=[mfilename('fullpath')];
    newbackup=strcat(outputPath,'/main_ga_backup.m');
    % newbackup=sprintf('%sbackup.m',FileNameAndLocation);
    currentfile=strcat(FileNameAndLocation, '.m');
    copyfile(currentfile,newbackup);
end















