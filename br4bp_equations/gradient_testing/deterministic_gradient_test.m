clc; close all; clear;

cd('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust')
addpath(genpath('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust'));

savename = ['br4bp_testing'];
scenario = 'deterministic br4bp testing';
saveOutput = true; % bool for saving the output or not, true or false
saveVideo = true;

% Setup for saving data
formatOut = 'yyyymmdd_HHMM.SS';
dateString = datestr(now,formatOut);
outputPath = strcat('./sims/',dateString,'_',savename);
%%

% init_fn = './init_traj_files/init_br4bp_x0_sb1_BLTii_2';
init_fn = './init_traj_files/init_br4bp_x0_sb1_BLTii_2p1';
% init_fn = './init_traj_files/init_br4bp_x0_sb1_BLTii_3';

run(init_fn);

%% 

traj0 = createStateStmHistory(simparams.x0, simparams);

[deltaV0, deltaVs_nom0] = calcDeltaV(simparams.x0, traj0.x_i_f, traj0.stm_i, simparams);








figure;
if simparams.perform_correction    
    plotMultiSegTraj(simparams.x0, traj0.x_t, traj0.t_s, simparams, tcm_idx0);
else
    plotMultiSegTraj(simparams.x0, traj0.x_t, traj0.t_s, simparams);
end
axis equal;



%% Fmincon call via output function
tic
mkdir(outputPath);

[x_opt,J,history,searchdir,exitflag,output] = runfmincon(simparams, outputPath);
toc

traj = createStateStmSttQdQHistory(x_opt, simparams);

figure
plotMultiSegTraj(x_opt, traj.x_t, traj.t_s, simparams);

%% 

if saveOutput
    saveallfigs(strcat('./',outputPath),0);
    save(strcat('./',outputPath,'/workspace.mat'));
end

% 
% 
% simparams.x0(3) = simparams.x0(3) + 1e-5;
% simparams.x0(5) = simparams.x0(5) + 1e-3;
% 
% 
% simparams.x0(end-16+3) = simparams.x0(end-16+3) + 1e-5;
% simparams.x0(end-16+6) = simparams.x0(end-16+6) + 1e-3;
% 
% 
% simparams.x0(end-5) = simparams.x0(end-5) + 1e-5;
% simparams.x0(end-2) = simparams.x0(end-2) + 1e-3;

% verifyGradients;
verifyGradients_circx0xf;

% % % x_opt(end) = 1e-16;


% Taking a really long time to propagate the seg in lunar orbit:
% removing the final segment
% simparams.n = 7;
% traj  = createStateStmSttHistory(simparams.x0(1:56), simparams);
% traj = createStateStmSttQdQHistory(simparams.x0(1:56), simparams);
% plotMultiSegTraj(simparams.x0(1:56), traj0.x_t, traj0.t_s, simparams);

traj = createStateStmSttQdQHistory(x_opt, simparams);
[deltaV, deltaVs_nom] = calcDeltaV(x_opt, traj.x_i_f, traj.stm_i, simparams);



[tcm_time, tcm_idx, min_tcm_dv, ~, ~, tcm_dv_each] = opt_multiple_tcm_wQ(x_opt, traj, deltaVs_nom, simparams);





solfig = figure;
% plotMultiSegTraj(x_opt, traj.x_t, traj.t_s, simparams, tcm_idx);
plotMultiSegTraj(x_opt, traj.x_t, traj.t_s, simparams);
title('optim_solution')
solfig.CurrentAxes.Title.Visible="off";

% % % 
% % % x = reshape(simparams.x0,simparams.m,simparams.n);
% % % % dx = sqrt(eps);
% % % dt = x(8,end);
% % % [~, stm] = stateStmProp(x(1:7,end), dt, simparams);
% % % 
% % % % Test the STT in the 8th segment
% % % for i = 1:7
% % %     xfd = x(1:7,end);
% % %     dx = 1e-5 * xfd(i);
% % %     
% % %     xfd(i) = xfd(i) + dx;
% % % 
% % %     [~, stm_fd] = stateStmProp(xfd, dt, simparams);
% % %     stt_fd(:,:,i) = (stm_fd - stm) ./dx;
% % % 
% % % 
% % % 
% % % 
% % % end




% % % % % for development purposes
% % % % simparams.tcm_nodes = [4, 5];
% % % % [cin, ceq, cinGrad, ceqGrad] = constraint_min_tcm(simparams.x0, simparams);



% Functions to modify for nsv
% calc_covariance_wQ_tcmdv_v3
% The rest of constraint_min_tcm when target_Pr_constraint_on = 1
