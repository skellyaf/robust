% 27 Apr 23: Summary of changes:
%       - Capability for additional nominal maneuvers. Current use cases are
%       for elliptical departures/flyby trajectories/a lap to clean up LV
%       dispersions.
%       - Mid-trajectory TCM target / covariance constraint not at the
%       target. Use case is to have a second nominal maneuver (TLI, for
%       example) where it is also desired that dispersions are minimized.
%       The TCM calculations are separated into two parts, before and after
%       the mid-course TCM target.
%       - Maneuver execution error. Implementing this required defining
%       covariance-modifying trajectory events, not just TCMs which were
%       the only events that modified the dispersion covariance previously.
%       This applies to both TCM maneuver execution error and nominal
%       maneuver execution error. Capability is also included that allows
%       them to occur at the same time.
%       - TCM improvement threshold. Rather than adding TCMs as long as the
%       total TCM delta V magnitude decreased, created a minimum value by
%       which it must decrease by (currently 3 sigma of the TCM execution
%       error).

% 1 Mar 23: Incorporating flexibility to CR3BP systems as well as 2BP.

% Analytical objective function and constraint gradients working

% Started new repo from multiseg_opt_manual_tcm on 19 jan 2023: purpose is
% to incorporate TCM execution error and multiple TCMs

clear;
clear global;
close all; 
clc;
format longg;
% cd('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust');
addpath(genpath('./'));

savename = 'cr3bp_heo_mlo_3dv_midTcmTgt_ndvErr_stdx0_robust';
saveOutput = true; % bool for saving the output or not, true or false


%% Create simulation parameters structure by running initialization script
% 2bp
% init_fn = './init_traj_files/initialize_simulation_parameters_45leo_to_45leo';
% init_fn = './init_traj_files/initialize_simulation_parameters_45leo_to_45leo_extra_coast';
% init_fn = './init_traj_files/initialize_simulation_parameters_leo28_to_geo0_targetEnd';
% init_fn = './init_traj_files/initialize_simulation_parameters_leo28_to_geo0';
% init_fn = './init_traj_files/initialize_simulation_parameters_leoCrit_to_leoSunSync';
% init_fn = './init_traj_files/initialize_simulation_parameters_45leo_to_45_20k';
% init_fn = './init_traj_files/initialize_simulation_parameters_45leoMoreCoast_to_45_20k';




% cr3bp
% init_fn = './init_traj_files/init_simparams_cr3bp_leo_to_llo';
% init_fn = './init_traj_files/init_simparams_cr3bp_geo_to_llo';
% init_fn = './init_traj_files/init_simparams_cr3bp_heo_to_mlo';
% init_fn = './init_traj_files/init_simparams_cr3bp_llo_to_nrho';

% cr3bp, 3 nominal maneuvers
% init_fn = './init_traj_files/init_simparams_cr3bp_heo_to_mlo_3dv';
init_fn = './init_traj_files/init_simparams_cr3bp_leo_to_mlo_3dv';
% init_fn = './init_traj_files/init_simparams_cr3bp_bigHeo_to_mlo_3dv';


run(init_fn);

%% ONLY If using a previous reference trajectory as the initial guess:
% load('.\init_traj_files\initial_guesses\tli_llo_deterministic_opt.mat')
% load('.\init_traj_files\initial_guesses\polar_llo_to_nrho_apolune.mat')

% simparams.x0 = x_opt;


%% Initialize transfer segments
% simparams = generateInitialXfer(simparams); % only for 2bp

%% View problem setup

% [~, x_i_f0, x_t0, stm_t0, t0, t_s0] = createStateStmHistory(simparams.x0, simparams);
[stm_i0, stt_i0, x_i_f0, x_t0, stm_t0, stt_t_i, t0, t_s0, stm_t_i0]  = createStateStmSttHistory(simparams.x0, simparams);

% Calculate total impulsive delta V for initial guess trajectory
[deltaV0, deltaVs_nom0] = calcDeltaV(simparams.x0,x_i_f0,simparams);
% [tcm_min, tcm_time, tcm_r, tcm_v, ~, ~, ~, ~, tcm_total_t] = tcmPair_rv(x_opt, t, stm_t, deltaVs_nom, simparams);


[tcm_time0, tcm_idx0, min_tcm_dv0, ~, ~, tcm_dv_each0] = opt_multiple_tcm(simparams.x0, t0, t_s0, stm_t0, simparams); % inputs: x, t, t_s, stm_t, stm_t_i, simparams
 

totalDV0 = deltaV0 + 3*min_tcm_dv0

figure;
if simparams.perform_correction    
    [~,tcm_idx0] = opt_multiple_tcm(simparams.x0, t0, t_s0, stm_t0, simparams);
    plotMultiSegTraj(simparams.x0, x_t0, t_s0, simparams, tcm_idx0);
else
    plotMultiSegTraj(simparams.x0, x_t0, t_s0, simparams);
end
axis equal;




%% Finite difference testing / gradient comparison

% verifyGradients;
% test_tcm_time_gradient;

%% TCM analysis development along initial guess trajectory
% multiple_tcm_min_r5;

%% Fmincon call via output function
tic
[x_opt,J,history,searchdir,exitflag,output] = runfmincon(simparams);
toc


%% View optimal solution

% Calculate optimal TCM time and delta V values
% [stm_i, x_i_f, x_t, stm_t, t, t_s] = createStateStmHistory(x_opt, simparams);
[stm_i, stt_i, x_i_f, x_t, stm_t, stt_t_i, t, t_s, stm_t_i]  = createStateStmSttHistory(x_opt, simparams);
% Calculate total impulsive delta V 
[deltaV, deltaVs_nom] = calcDeltaV(x_opt, x_i_f, simparams);
% [tcm_min, tcm_time, tcm_r, tcm_v, ~, ~, ~, ~, tcm_total_t] = tcmPair_rv(x_opt, t, stm_t, deltaVs_nom, simparams);
[tcm_time,tcm_idx,min_tcm_dv] = opt_multiple_tcm(x_opt, t, t_s, stm_t, simparams);
totalDV = deltaV + 3*min_tcm_dv



solfig = figure;
plotMultiSegTraj(x_opt, x_t, t_s, simparams, tcm_idx);
title('optim_solution')
solfig.CurrentAxes.Title.Visible="off";


% v2 = x_opt(11:13);
% [~,~,i] = rv2orbel(r2,v2,simparams.mu);
% xfer_inclination = i*180/pi




%% Save outputs

% if output.firstorderopt < 4 && exitflag ~= -2 && saveOutput
if exitflag ~= -2 && saveOutput


    % Saving movie of history
    formatOut = 'yyyymmdd_HHMM.SS';
    dateString = datestr(now,formatOut);
    outputPath = strcat('./sims/',dateString,'_',savename);
    mkdir(outputPath);
    outputPathName = strcat('./',outputPath,'/video.avi');

    % Save workspace    
    save(strcat('./',outputPath,'/workspace.mat'));
    saveallfigs(strcat('./',outputPath),0)
    
    % save main file backup
    FileNameAndLocation=[mfilename('fullpath')];
    newbackup=strcat(outputPath,'/main_optim_backup.m');
    % newbackup=sprintf('%sbackup.m',FileNameAndLocation);
    currentfile=strcat(FileNameAndLocation, '.m');
    copyfile(currentfile,newbackup);
    
    % save constraint function backup
    newbackup=strcat(outputPath,'/constraint_min_tcm_backup.m');
    copyfile('./constraint_min_tcm.m',newbackup);
    
    % save objective function backup
    newbackup=strcat(outputPath,'/obj_min_tcm_backup.m');
    copyfile('./obj_min_tcm.m',newbackup);
    
    % save tcm optimization function backup
    newbackup=strcat(outputPath,'/tcmPair_rv_backup.m');
    copyfile('./Analysis/tcmPair_rv.m',newbackup);
    
    % save simulation parameters backup
%     newbackup=strcat('./',outputPath,'/initialize_simulation_parameters_backup.m');
    newbackup=strcat(outputPath,'/initialize_simparams_backup.m');
    copyfile(strcat(init_fn,'.m'),newbackup);



    optimHistoryMovie(history,outputPathName, simparams);


end




