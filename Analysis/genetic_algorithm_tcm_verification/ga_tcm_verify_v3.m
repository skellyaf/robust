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


savename = 'nri_flex0_robustverify';

%% x0

% [stm_i0, stt_i0, x_i_f0, x_t0, stm_t0, stt_t_i, t0, t_s0, stm_t_i0]  = createStateStmSttHistory(simparams.x0, simparams);
% [tcm_time0, tcm_idx0, min_tcm_dv0, ~, ~, tcm_dv_each0] = opt_multiple_tcm(simparams.x0, t0, t_s0, stm_t0, simparams); % inputs: x, t, t_s, stm_t, stm_t_i, simparams



% commented this out when loading a workspace directly below
% traj  = createStateStmSttQdQHistory(x, simparams);
% [deltaV, deltaVs_nom] = calcDeltaV(x, traj.x_i_f, traj.stm_i, simparams);
% 
% [tcm_time, tcm_idx, min_tcm_dv, ~, ~, tcm_dv_each] = opt_multiple_tcm_wQ(x, traj, deltaVs_nom, simparams); % inputs: x, t, t_s, stm_t, stm_t_i, simparams

% load('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\20240228_1754.43_3dv_nri_meddx0_flybynotcorrected_robust\workspace.mat')


% load('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\20240229_0837.29_3dv_nri_meddx0_flybynotcorrected_det\workspace.mat')

% load('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\20240222_1727.03_2dv_leo_llo_det\workspace.mat')

% det opt nri flex0
load('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\20240228_1252.03_3dv_nri_smalldx0_flybynotcorrected_det_flex0\workspace.mat')

% init_fn = 'init_simpar_nri_3dv_TCMsNodes_flybyNotCorrected_flex0';
% init_fn = 'init_simpar_nri_3dv_TCMsNodes_flybyNotCorrected_startflex0';
% run(init_fn);


% robust flex0
% load('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\20240316_1624.04_3dv_nri_meddx0_flybynotcorrected_robust_flex0\workspace.mat')
% load('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\from_aries\20240317_1711.31_3dv_nri_flybynotcorrected_robust_flex0\workspace.mat')
% load('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\20240316_2040.59_3dv_nri_meddx0_flybynotcorrected_robust_flex0\workspace.mat')
% load('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\from_aries\20240317_1711.54_3dv_nri_flybynotcorrected_robust_flex0\workspace.mat')
% load('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\from_aries\20240317_2021.10_3dv_nri_flybynotcorrected_robust_flex0\workspace.mat')


% close all;

x = reshape(x_opt, simparams.m, simparams.n);




%% IF CHANGING THE STOCHASTIC PARAMS

simparams.P_max_r = 1 / ndDist2km; % km converted to ND dist
% 
% % Initial uncertainty
% Very Small
% simparams.sig_pos = 10 / 1e3 / ndDist2km; % Position +/- 10 m in all 3 direction
% simparams.sig_vel = 10 / 1e6 / ndDist2km * ndTime2sec; % Velocity +/- 1 cm/s in all 3 directions

% Small
% simparams.sig_pos = 100 / 1e3 / ndDist2km; % Position +/- 100 m in all 3 direction
% simparams.sig_vel = .1 / 1e3 / ndDist2km * ndTime2sec; % Velocity +/- 10 cm/s in all 3 directions

% Medium
simparams.sig_pos = 1 / ndDist2km; % Position +/- 1 km in all 3 direction converted to ND dist
simparams.sig_vel = 1 / 1e3 / ndDist2km * ndTime2sec; % Velocity +/- 1 m/s in all 3 directions converted to ND dist / ND time

simparams.P_initial = diag([simparams.sig_pos^2 simparams.sig_pos^2 simparams.sig_pos^2 simparams.sig_vel^2 simparams.sig_vel^2 simparams.sig_vel^2]);

% TCM execution error
% simparams.sig_tcm_error = .01 / 1e3 / ndDist2km * ndTime2sec; % Velocity 1 sigma = 1 cm/s
simparams.sig_tcm_error = .1 / 1e3 / ndDist2km * ndTime2sec; % Velocity 1 sigma = 10 cm/s
simparams.R = diag([simparams.sig_tcm_error, simparams.sig_tcm_error, simparams.sig_tcm_error]).^2;

% Nominal maneuver execution error
% simparams.sig_dv_error = 1 / 1e3 / ndDist2km * ndTime2sec; % Velocity 1 sigma = 1 m/s
simparams.sig_dv_error = 10 / 1e3 / ndDist2km * ndTime2sec; % Velocity 1 sigma = 10 m/s
simparams.R_dv = diag([simparams.sig_dv_error, simparams.sig_dv_error, simparams.sig_dv_error]).^2;

simparams.add_tcm_improvement_threshold = sqrt(trace(simparams.R)) * 3;


% simparams.Qt = 1e-8 * eye(3) * ndTime2sec^3 / ndDist2m^2; % m^2 / sec^3 converted to ND dist ^ 2 / ND time ^ 3
simparams.Qt = 1e-6 * eye(3) * ndTime2sec^3 / ndDist2m^2; % m^2 / sec^3 converted to ND dist ^ 2 / ND time ^ 3
% simparams.Qt = (.3e-3)^2 * eye(3) * ndTime2sec^3 / ndDist2m^2; % m^2 / sec^3 converted to ND dist ^ 2 / ND time ^ 3
% simparams.Qt = 1e-5 * eye(3) * ndTime2sec^3 / ndDist2m^2; % m^2 / sec^3 converted to ND dist ^ 2 / ND time ^ 3



% %%
traj = createStateStmSttQdQHistory(x(:), simparams);
[deltaV, deltaVs_nom] = calcDeltaV(x, traj.x_i_f, traj.stm_i, simparams);
% % IF WANT TO TEST THE RESULT OF THE ALGORITHM ALONG THE NOMINAL TRAJ --
% % UNCOMMENT IF TESTING THE NLP TCMs
[tcm_time, tcm_idx, min_tcm_dv, ~, ~, tcm_dv_each] = opt_multiple_tcm_wQ_multiPart(x(:), traj, deltaVs_nom, simparams); % inputs: x, t, t_s, stm_t, stm_t_i, simparams
[Q_k_km1, dQ_k_km1_dxi, dQ_k_km1_ddti] = calc_Q_events(traj, x_opt, tcm_time, simparams);

[P_target, min_tcm_dv, tcm_dv, P_i_minus, P_i_plus] = calc_covariance_wQ_tcmdv_v3(x, traj, tcm_time, 1, deltaVs_nom, simparams.P_initial, Q_k_km1, simparams);


min_tcm_dv*ndVel2kms*3000




%%


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


    %%%%%%%%%%% REMOVE AFTER TEST %%%%%%%%%%%%%%%%
%     num_tcms_p2 = num_tcms_p2 - 1;


end

% num_tcms_p1 = num_tcms_p1 + 1;
% num_tcms_p2 = num_tcms_p2+1;


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

% options = optimoptions(@ga,...
%             		'ConstraintTolerance',1e-12,...
%             		'FunctionTolerance',1e-12,...
%             		'MaxTime',36000,...
%            		    'MaxGenerations',500 ,...
%            		    'Display','iter',...
%             		'UseParallel',true, ...
%                     'InitialPopulationMatrix',tcm_idx', ...
%                     'PopulationSize',1, ...
%                     'MaxStallGenerations',50);

options = optimoptions(@ga,...
            		'ConstraintTolerance',1e-12,...
            		'FunctionTolerance',1e-12,...
            		'MaxTime',36000,...
           		    'MaxGenerations',500 ,...
           		    'Display','iter',...
            		'UseParallel',true, ...
                    'PopulationSize',100, ...
                    'MaxStallGenerations',50);



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
%     num_tcms = length(tcm_time);
    num_tcms = num_tcms_p1 + num_tcms_p2;
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















