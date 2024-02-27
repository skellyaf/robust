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
addpath(genpath('./'));

savename = ['3dv_eed_llo_det'];
scenario = '3dv LEO to EED to LLO deterministic';
saveOutput = true; % bool for saving the output or not, true or false
saveVideo = true;

% Setup for saving data
formatOut = 'yyyymmdd_HHMM.SS';
dateString = datestr(now,formatOut);
outputPath = strcat('./sims/',dateString,'_',savename);
    


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
% init_fn = './init_traj_files/init_simparams_cr3bp_leo_to_llo_loadDetOpt'; %%%%%%%%
% init_fn = './init_traj_files/init_simparams_cr3bp_geo_to_llo';
% init_fn = './init_traj_files/init_simparams_cr3bp_heo_to_mlo';
% init_fn = './init_traj_files/init_simparams_cr3bp_llo_to_nrho';

% cr3bp, 3 nominal maneuvers
% init_fn = './init_traj_files/init_simparams_cr3bp_heo_to_mlo_3dv';
% init_fn = './init_traj_files/init_simparams_cr3bp_leo_to_mlo_3dv';

% init_fn = './init_traj_files/init_simparams_cr3bp_leo_lloflyby_nri_3dv';
% init_fn = './init_traj_files/init_simparams_cr3bp_leoinclined_lloflyby_nri_3dv';

% init_fn = './init_traj_files/init_simparams_cr3bp_leoinclined_lloflyby_dro_3dv';
% init_fn = './init_traj_files/init_simparams_cr3bp_leoinclined_lloflyby_dro2_3dv';
% init_fn = './init_traj_files/init_simparams_cr3bp_leoinclined_lloflyby_dro3_3dv';
% init_fn = './init_traj_files/init_simparams_cr3bp_leoinclined_lloflyby_dro4_3dv';

% init_fn = './init_traj_files/init_simparams_cr3bp_leoinclined_lloflyby_dro3_x0detOpt_3dv';


% 3 dv NRI, TCMs at nodes
% init_fn = './init_traj_files/init_simparams_cr3bp_leoinclined_lloflyby_nri_3dv_TCMsAtNodes';
% init_fn = './init_traj_files/init_simpar_leoinclined_lloflyby_nri_3dv_TCMsNodes_moreQ';
% init_fn = './init_traj_files/init_simpar_nri_3dv_TCMsNodes_flybyNotCorrected';

% EED to LEO, TCMs at nodes
init_fn = './init_traj_files/init_simparams_cr3bp_leo_to_mlo_3dv_tcmAtNodes'; %%%%%%%%%%%%%%%%
% init_fn = './init_traj_files/init_simparams_cr3bp_leo_to_mlo_3dv_tcmAtNodes_loadSimparams';


% init_fn = './init_traj_files/init_simparams_cr3bp_bigHeo_to_mlo_3dv';
% init_fn = './init_traj_files/init_simparams_cr3bp_flex0_leo_lloflyby_nri_3dv';


% Rendezvous, nrho
% init_fn = './init_traj_files/init_simparams_cr3bp_nrho_rdvz_2dv';

% Rendezvous, nrho, make the TCMs occur at segment intersections
% init_fn = './init_traj_files/init_sim_nrho_rdvz_2dv_tcmAtSegIntersection'; %%%%%%%


run(init_fn);




%%%%% DEBUGGING

% load('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\20231010_1041.47_nri_3dv_TcmAtNodes_robust_WORKED\workspace.mat')
% simparams.sig_dv_error = 10 / 1e3 / ndDist2km * ndTime2sec; % Velocity 1 sigma = 10 m/s
% simparams.R_dv = diag([simparams.sig_dv_error, simparams.sig_dv_error, simparams.sig_dv_error]).^2;
% simparams.constrain_flyby_radius = 0;
% formatOut = 'yyyymmdd_HHMM.SS';
% dateString = datestr(now,formatOut);
% outputPath = strcat('./sims/',dateString,'_',savename);
% mkdir(outputPath);

%%%%%%%%



%% ONLY If using a previous reference trajectory as the initial guess:
% load('.\init_traj_files\initial_guesses\tli_llo_deterministic_opt.mat')
% load('.\init_traj_files\initial_guesses\polar_llo_to_nrho_apolune.mat')

% load('eed_leo_planar_13day.mat');
% load('nri_det_opt.mat');
% load('nri_det_opt_20seg.mat');
% load('nri_planar_det_opt.mat');
% load('leo_plf_dro3_detOpt.mat')

% load('nrho_rdvz_det_opt_12hrlead.mat');
% load('nrho_rdvz_gt1orbit.mat');
% 
% simparams.x0 = x_opt(:);


%% Initialize transfer segments
% simparams = generateInitialXfer(simparams); % only for 2bp

%% View problem setup


% [~, x_i_f0, x_t0, stm_t0, t0, t_s0] = createStateStmHistory(simparams.x0, simparams);
% [stm_i0, stt_i0, x_i_f0, x_t0, stm_t0, stt_t_i, t0, t_s0, stm_t_i0]  = createStateStmSttHistory(simparams.x0, simparams);
traj0  = createStateStmSttQdQHistory(simparams.x0, simparams);

if isfield(simparams,'rdvz_flag')
    if simparams.rdvz_flag == 1
        x0 = reshape(simparams.x0,simparams.m,simparams.n);
        totalTime = sum(x0(7,:));
        simparams.x_target = stateProp(simparams.x0_target, totalTime, simparams);
    end
end

% Calculate total impulsive delta V for initial guess trajectory
[deltaV0, deltaVs_nom0] = calcDeltaV(simparams.x0, traj0.x_i_f, traj0.stm_i, simparams);


% [tcm_time0, tcm_idx0, min_tcm_dv0, ~, ~, tcm_dv_each0] = opt_multiple_tcm(simparams.x0, deltaVs_nom0, t0, t_s0, stm_t0, stm_t_i0, simparams); % inputs: x, t, t_s, stm_t, stm_t_i, simparams
% [tcm_time0, tcm_idx0, min_tcm_dv0, ~, ~, tcm_dv_each0] = opt_multiple_tcm_wQ(simparams.x0, traj0, deltaVs_nom0, simparams); % inputs: x, t, t_s, stm_t, stm_t_i, simparams

if ~simparams.perform_correction
    [tcm_time0, tcm_idx0, min_tcm_dv0, ~, ~, tcm_dv_each0] = opt_multiple_tcm_wQ(simparams.x0, traj0, deltaVs_nom0, simparams); % inputs: x, t, t_s, stm_t, stm_t_i, simparams
else
    tcm_time0 = zeros(1, length(simparams.tcm_nodes));
    
    x0 = reshape(simparams.x0, simparams.m, simparams.n);
    for i = 1:length(simparams.tcm_nodes)                
        tcm_time0(i) = sum(x0(7,1:simparams.tcm_nodes(i)-1));
    end
    
    
    for i = 1:length(tcm_time0)
        tcm_idx0(i) = find(traj0.t == tcm_time0(i))';
    end

    

end


[Q_k_km10, dQ_k_km1_dxi0, dQ_k_km1_ddti0] = calc_Q_events(traj0, simparams.x0, tcm_time0, simparams);

[P_target0, min_tcm_dv0, tcm_dv_each0, P_i_minus0, P_i_plus0] = calc_covariance_wQ_tcmdv_v3(simparams.x0, traj0, tcm_time0, 1, deltaVs_nom0, simparams.P_initial, Q_k_km10, simparams);


totalDV0 = deltaV0 + 3*min_tcm_dv0

disp('Initial guess with optimal corrections')
deltaV0*ndVel2kms*1000, 3*min_tcm_dv0*ndVel2kms*1000, totalDV0*ndVel2kms*1000

figure;
if simparams.perform_correction    
    plotMultiSegTraj(simparams.x0, traj0.x_t, traj0.t_s, simparams, tcm_idx0);
else
    plotMultiSegTraj(simparams.x0, traj0.x_t, traj0.t_s, simparams);
end
axis equal;




%% Finite difference testing / gradient comparison

% verifyGradients;
% test_tcm_time_gradient;

%% TCM analysis development along initial guess trajectory
% multiple_tcm_min_r5;

%% Fmincon call via output function
tic
mkdir(outputPath);

[x_opt,J,history,searchdir,exitflag,output] = runfmincon(simparams, outputPath);
toc


%% View optimal solution

% Calculate optimal TCM time and delta V values
% [stm_i, stt_i, x_i_f, x_t, stm_t, stt_t_i, t, t_s, stm_t_i]  = createStateStmSttHistory(x_opt, simparams);
traj = createStateStmSttQdQHistory(x_opt, simparams);


if isfield(simparams,'rdvz_flag')
    if simparams.rdvz_flag == 1
        x_opt = reshape(x_opt,simparams.m,simparams.n);
        totalTime = sum(x_opt(7,:));
        simparams.x_target = stateProp(simparams.x0_target, totalTime, simparams);
    end
end

% Calculate total impulsive delta V 
[deltaV, deltaVs_nom] = calcDeltaV(x_opt, traj.x_i_f, traj.stm_i, simparams);
% [tcm_time,tcm_idx,min_tcm_dv,~,~,tcm_dv_each] = opt_multiple_tcm(x_opt, deltaVs_nom, t, t_s, stm_t, stm_t_i, simparams);
if ~simparams.perform_correction
    [tcm_time, tcm_idx, min_tcm_dv, ~, ~, tcm_dv_each] = opt_multiple_tcm_wQ(x_opt, traj, deltaVs_nom, simparams); % inputs: x, t, t_s, stm_t, stm_t_i, simparams
else
    tcm_time = zeros(1, length(simparams.tcm_nodes));
    
    x_opt = reshape(x_opt, simparams.m, simparams.n);
    for i = 1:length(simparams.tcm_nodes)                
        tcm_time(i) = sum(x_opt(7,1:simparams.tcm_nodes(i)-1));
    end
    
    
    for i = 1:length(tcm_time)
        tcm_idx(i) = find(traj.t == tcm_time(i))';
    end
end


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

[event_times, event_indicator] = define_events_v2(x_opt(:), traj.t, tcm_time, simparams);






%%%% views of the flex0 TLI traj

solfig2 = figure;
plotMultiSegTraj(x_opt, traj.x_t, traj.t_s, simparams, tcm_idx);
% plotMultiSegTraj(x_opt, x_t, t_s, simparams);
axis equal;
title('optim_solution_view2')
solfig2.CurrentAxes.Title.Visible="off";


xlim([.9 1.025])
zlim([-.01 .05])
ylim([-.025 .15])



% Comparison with 'LAR' TCM method -- for NRHO example for now
% x_opt = reshape(x_opt,simparams.m,simparams.n);
% t_flyby = sum(x_opt(7,1:simparams.P_constrained_nodes(1)));
% t_nri = sum(x_opt(7,1:simparams.P_constrained_nodes(2)));
% 
% 
% 
% tcm_idx_lar = [1500 2200 tcm_idx(5), 6000, tcm_idx(end)];
% tcm_time_lar = traj.t(tcm_idx_lar)';
% 
% Q_k_km1 = calc_Q_events(traj, x_opt, tcm_time_lar, simparams);
% 
% [P_lar, min_tcm_dv_total_lar, tcm_dv_each_lar, P_i_minus_lar] = calc_covariance_wQ_tcmdv_v3(x_opt(:), traj, tcm_time_lar, 1, deltaVs_nom, simparams.P_initial, Q_k_km1, simparams);
% larfig = figure;
% plotMultiSegTraj(x_opt, traj.x_t, traj.t_s, simparams, tcm_idx_lar);
% title('opt_sol_w_lar_tcms')
% larfig.CurrentAxes.Title.Visible="off";
% 
% 
% disp('Robust with LAR corrections')
% 3*min_tcm_dv_total_lar*ndVel2kms*1000

%% Save iteration history 
iterhist = figure;
plotIterationHistory(history.x,simparams)
title('iteration_history')
iterhist.CurrentAxes.Title.Visible="off";


xlim([0.98 1.02])
ylim([0 .05])
zlim([-.1375, -.075])

%% Save outputs

% if output.firstorderopt < 4 && exitflag ~= -2 && saveOutput
% if exitflag ~= -2 && saveOutput
if saveOutput


    
 
    % Save workspace    
    saveallfigs(strcat('./',outputPath),0);
    save(strcat('./',outputPath,'/workspace.mat'));
    
    
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
    newbackup=strcat(outputPath,'/opt_multiple_tcm_backup.m');
    copyfile('./Analysis/multiple_tcm_devt/opt_multiple_tcm.m',newbackup);

    % Analysis/multiple_tcm_devt/min_dv_tcm_meets_dispersion_constraint
    newbackup=strcat(outputPath,'/min_dv_tcm_meets_dispersion_constraint_backup.m');
    copyfile('./Analysis/multiple_tcm_devt/min_dv_tcm_meets_dispersion_constraint.m',newbackup);

    % Analysis/multiple_tcm_devt/optimize_tcm_guess
    newbackup=strcat(outputPath,'/optimize_tcm_guess_backup.m');
    copyfile('./Analysis/multiple_tcm_devt/optimize_tcm_guess.m',newbackup);

    % Analysis/multiple_tcm_devt/tcm_index_gradient_vector_search
    newbackup=strcat(outputPath,'/tcm_index_gradient_vector_search_backup.m');
    copyfile('./Analysis/multiple_tcm_devt/tcm_index_gradient_vector_search.m',newbackup);

    % Analysis/multiple_tcm_devt/random_unit_search
    newbackup=strcat(outputPath,'/random_unit_search_backup.m');
    copyfile('./Analysis/multiple_tcm_devt/random_unit_search.m',newbackup);

    % Analysis/calc_covariance_tcmdv
    newbackup=strcat(outputPath,'/calc_covariance_tcmdv_backup.m');
    copyfile('./Analysis/calc_covariance_tcmdv.m',newbackup);

    
    % save simulation parameters backup
    newbackup=strcat(outputPath,'/initialize_simparams_backup.m');
    copyfile(strcat(init_fn,'.m'),newbackup);



    if saveVideo
        vidOutputPathName = strcat('./',outputPath,'/video.avi');
        optimHistoryMovie(history,vidOutputPathName, simparams);
    end


end



% savefig(iterhist, strcat(outputPath,'/iteration_history.fig'));

% exportgraphics(iterhist,'iterhist.pdf','ContentType','image')



%% Save to results summary file

if 1
    save_to_excel;
%     save_to_excel_nrho_rdvz;
    
    

end

%% debug
[traj0]  = createStateStmSttHistory(x, simparams);
% [tcm_time0, tcm_idx0, min_tcm_dv0, ~, ~, tcm_dv_each0] = opt_multiple_tcm(x, t0, t_s0, stm_t0, simparams); % inputs: x, t, t_s, stm_t, stm_t_i, simparams
figure
plotMultiSegTraj(x, traj0.x_t, traj0.t_s, simparams);




[traj0]  = createStateStmSttHistory(x_opt, simparams);
% [tcm_time0, tcm_idx0, min_tcm_dv0, ~, ~, tcm_dv_each0] = opt_multiple_tcm(x, t0, t_s0, stm_t0, simparams); % inputs: x, t, t_s, stm_t, stm_t_i, simparams
figure
plotMultiSegTraj(x_opt, traj0.x_t, traj0.t_s, simparams);



%% investigating natural motion of 


% [~, x_nofinalDV_t] = stateProp(x_opt(1:6,end-1),50/ndTime2days,simparams);
% 
% figure
% axis equal
% grid on
% hold on
% plot3(x_nofinalDV_t(:,1),x_nofinalDV_t(:,2),x_nofinalDV_t(:,3),'LineWidth',2)

%% investigating apse constraint
% v_I = zeros(3,length(t));
% 
% rdotv = zeros(1,length(t));
% r_m_sc_mag = rdotv;
% 
% for i = 1:length(t)
% 
%     omega_SI_nd = [0; 0; 1];
%     r_m = [1-mu; 0; 0];
% 
%     % Get inertial velocity
%     
%     v_S = x_t(i,4:6)';
%     r_S = x_t(i,1:3)';
%     v_I(:,i) = v_S + cross(omega_SI_nd, r_S);
% 
% 
%     r_m_sc = r_S - r_m;
%     r_m_sc_mag(i) = vecnorm(r_m_sc);
%     i_r_m_sc = r_m_sc/vecnorm(r_m_sc);
% 
%     i_v_I = v_I(:,i)/vecnorm(v_I(:,i));
% 
%     rdotv(i) = i_r_m_sc' * i_v_I;
% 
%     
% end
% 
% 
% figure
% plot(rdotv)
% hold on;
% yline(0)
% plot(r_m_sc_mag)


% 
% test_dot = abs(rdotv(4100:4300));
% [minv,minidx] = min(test_dot);
% 
% rdotv(4100+minidx-1)
% 
% 
% 
% plot3(x_t(4100+minidx-1,1),x_t(4100+minidx-1,2),x_t(4100+minidx-1,3),'.','MarkerSize',25)


beep