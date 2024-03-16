%% dev
clear all;
clc;


%% Load previous workspace for solution

% load('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\20240229_0837.29_3dv_nri_meddx0_flybynotcorrected_det\workspace.mat')
% load('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\20240222_1727.03_2dv_leo_llo_det\workspace.mat')
% load('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\from_aries\20231016_1050.20_eed_llo_3dv_TcmAtNodes_det\workspace.mat')
% load('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\20240229_0837.29_3dv_nri_meddx0_flybynotcorrected_det\workspace.mat')
% load('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\20240228_1754.43_3dv_nri_meddx0_flybynotcorrected_robust\workspace.mat')
load('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\20240314_1510.58_3dv_eed_llo_det\workspace.mat')



savename = ['rerun_new_tcm_opt_leo_eed_llo_det'];

close all;


%% Random variables / initial uncertainty
% Units km, km/hr, km/hr^2

% simparams.P_max_r = 100 / ndDist2km; % km converted to ND dist
simparams.P_max_r = 1 / ndDist2km; % km converted to ND dist

% Initial uncertainty
% zero
% simparams.sig_pos = 1e-12;
% simparams.sig_vel = 1e-12;

% Small
simparams.sig_pos = 10 / 1e3 / ndDist2km; % Position +/- 10 m in all 3 direction
simparams.sig_vel = 10 / 1e6 / ndDist2km * ndTime2sec; % Velocity +/- 1 cm/s in all 3 directions

% Medium
% simparams.sig_pos = 1 / ndDist2km; % Position +/- 1 km in all 3 direction converted to ND dist
% simparams.sig_vel = 1 / 1e3 / ndDist2km * ndTime2sec; % Velocity +/- 1 m/s in all 3 directions converted to ND dist / ND time

% Large
% simparams.sig_pos = 10 / ndDist2km; % Position +/- 10 km in all 3 direction converted to ND dist
% simparams.sig_vel = 10 / 1e3 / ndDist2km * ndTime2sec; % Velocity +/- 10 m/s in all 3 directions converted to ND dist / ND time
% simparams.sig_vel = 10 / 1e5 / ndDist2km * ndTime2sec; % Velocity +/- 10 cm/s in all 3 directions converted to ND dist / ND time


simparams.P_initial = diag([simparams.sig_pos^2 simparams.sig_pos^2 simparams.sig_pos^2 simparams.sig_vel^2 simparams.sig_vel^2 simparams.sig_vel^2]);

% TCM execution error
% simparams.sig_tcm_error = 1e-12; % Velocity 1 sigma = nearly 0 cm/s
% simparams.sig_tcm_error = .01 / 1e3 / ndDist2km * ndTime2sec; % Velocity 1 sigma = 1 cm/s
simparams.sig_tcm_error = .1 / 1e3 / ndDist2km * ndTime2sec; % Velocity 1 sigma = 10 cm/s

simparams.R = diag([simparams.sig_tcm_error, simparams.sig_tcm_error, simparams.sig_tcm_error]).^2;

% Nominal maneuver execution error
% simparams.sig_dv_error = 1e-12; % Velocity 1 sigma = nearly 0 cm/s
% simparams.sig_dv_error = .1 / 1e3 / ndDist2km * ndTime2sec; % Velocity 1 sigma = 10 cm/s
simparams.sig_dv_error = 10 / 1e3 / ndDist2km * ndTime2sec; % Velocity 1 sigma = 10 m/s


simparams.R_dv = diag([simparams.sig_dv_error, simparams.sig_dv_error, simparams.sig_dv_error]).^2;


simparams.add_tcm_improvement_threshold = sqrt(trace(simparams.R)) * 3;

% simparams.R = diag([0 0 0]);

% simparams.Qt = sqrt(4.8e-7^2 / 3) * eye(3) * (ndTime2sec^3/ndDist2km^2) * .000001; % the value used for dev/testing
% simparams.Qt = sqrt(4.8e-7^2 / 3) * eye(3) * (ndTime2sec^3/ndDist2km^2) * .00001;
% simparams.Qt = 4.8e-7 * eye(3) * ndTime2sec^3 / ndDist2m^2 * 1;
simparams.Qt = 1e-6 * eye(3);
% simparams.Qt = 1e-8 * eye(3);


%% Find new optimal TCM solution
traj = createStateStmSttQdQHistory(x_opt, simparams);

[tcm_time, tcm_idx, min_tcm_dv, ~, ~, tcm_dv_each] = opt_multiple_tcm_wQ_multiPart(x_opt, traj, deltaVs_nom, simparams)





[Q_k_km1, dQ_k_km1_dxi, dQ_k_km1_ddti] = calc_Q_events(traj, x_opt, tcm_time, simparams);

[P_target, min_tcm_dv, tcm_dv, P_i_minus, P_i_plus] = calc_covariance_wQ_tcmdv_v3(x_opt, traj, tcm_time, 1, deltaVs_nom, simparams.P_initial, Q_k_km1, simparams);


totalDV = deltaV + 3*min_tcm_dv




disp('Optimal')
deltaVmps = deltaV*ndVel2kms*1000
tcmtotalmps = 3*min_tcm_dv*ndVel2kms*1000
totalcostmps = totalDV*ndVel2kms*1000




solfig = figure;
plotMultiSegTraj(x_opt, traj.x_t, traj.t_s, simparams, tcm_idx);
% plotMultiSegTraj(x_opt, x_t, t_s, simparams);
axis equal;
title('optim_solution')
solfig.CurrentAxes.Title.Visible="off";


%% Plot



%% Save



saveOutput = true; % bool for saving the output or not, true or false

% Setup for saving data
formatOut = 'yyyymmdd_HHMM.SS';
dateString = datestr(now,formatOut);
outputPath = strcat('./sims/',dateString,'_',savename);

mkdir(outputPath);
    
if saveOutput


    
 
    % Save workspace    
    saveallfigs(strcat('./',outputPath),0);
    save(strcat('./',outputPath,'/workspace.mat'));
    
    
    % save main file backup
    FileNameAndLocation=[mfilename('fullpath')];
    newbackup=strcat(outputPath,'/main_backup.m');
    % newbackup=sprintf('%sbackup.m',FileNameAndLocation);
    currentfile=strcat(FileNameAndLocation, '.m');
    copyfile(currentfile,newbackup);
    
    % Analysis/multiple_tcm_devt/min_dv_tcm_meets_dispersion_constraint
    newbackup=strcat(outputPath,'/min_dv_tcm_meets_dispersion_constraint_backup.m');
    copyfile('./Analysis/multiple_tcm_devt/min_dv_tcm_meets_dispersion_constraint_wQ_v3.m',newbackup);

    % Analysis/multiple_tcm_devt/optimize_tcm_guess
    newbackup=strcat(outputPath,'/optimize_tcm_guess_backup.m');
    copyfile('./Analysis/multiple_tcm_devt/optimize_tcm_guess_wQ_entireTraj.m',newbackup);

    % Analysis/multiple_tcm_devt/tcm_index_gradient_vector_search
    newbackup=strcat(outputPath,'/tcm_index_gradient_vector_search_backup.m');
    copyfile('./Analysis/multiple_tcm_devt/tcm_index_gradient_vector_search_wQ_entireTraj.m',newbackup);

    % Analysis/calc_covariance_tcmdv
    newbackup=strcat(outputPath,'/calc_covariance_tcmdv_backup.m');
    copyfile('./Analysis/calc_covariance_wQ_tcmdv_v3.m',newbackup);


end