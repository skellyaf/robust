clear;
clear global;
close all; 
clc;
format longg;
addpath(genpath('../../'));


% load('C:\Users\skell\OneDrive - USU\Documents\code_repos\multiseg_opt_manual_tcm\sims\20221103_1138.29_leo45_to_leoSunSync_testing_deterministic_incChangeAllDV1_targetXFN_freeTCMloc\workspace.mat')
% load('C:\Users\skell\OneDrive - USU\Documents\code_repos\multiseg_opt_manual_tcm\sims\20221101_1623.28_leo45_to_leoSunSync_robust_10kmR_incChangeAllDV1_J\workspace.mat')
load('C:\Users\skell\OneDrive - USU\Documents\code_repos\multiseg_opt_manual_tcm\sims\20221101_1333.01_leo28_to_geo0_J\workspace.mat')
simparams.sig_pos = 100;
simparams.sig_vel = .1;
simparams.P_initial = diag([simparams.sig_pos^2 simparams.sig_pos^2 simparams.sig_pos^2 simparams.sig_vel^2 simparams.sig_vel^2 simparams.sig_vel^2]);


%% Monte carlo simulation of trajectory x_opt

x_opt = reshape(x_opt,simparams.m, simparams.n);
% Simulate trajectory and determine miss distance at the target




% Propagate unperturbed initial trajectory, calculate dv1 and dv2

[stm_i, x_i_f, x_t, stm_t, t] = createStateStmHistory(x_opt, simparams);
[~, deltaVs_nom] = calcDeltaV(x_opt,x_i_f,simparams);

maneuver_times = [];
for i = 1:length(simparams.maneuverSegments)
    maneuver_times(i) = sum(x_opt(7,1:simparams.maneuverSegments(i)-1));
end

% Determine target states and times (maneuver/final trajectory...I suppose
% tcm_target is the only one that matters)
final_target = simparams.x_target;
final_target_time = sum(x_opt(7,:));

% simparams.target_final_maneuver = 1; %%%% CHANGED TO TARGET FINAL MANEUVER

if simparams.target_final_maneuver
    tcm_target = x_opt(1:6,simparams.maneuverSegments(end)) - [zeros(3,1); deltaVs_nom(:,2)];
    tcm_target_time = sum(x_opt(7,1:simparams.maneuverSegments(end)-1));
else
    tcm_target = final_target;
    tcm_target_time = final_target_time;
end

% Get the tcm execution time
[tcm_min,tcm_time] = tcmPair_rv(x_opt, t, stm_t, simparams);

%%%%%% CHANGED TO ALIGN TCM WITH FIRST DELTA V
% tcm_time = maneuver_times(1);


figure
plot3(x_t(:,1),x_t(:,2),x_t(:,3),'LineWidth',3,'Color','Black')
axis equal
hold on
grid on
plot3(tcm_target(1),tcm_target(2),tcm_target(3),'.','markersize',25)

tcm_norm_sav = [];

%% Start random variable part
for j = 1:200

    
    % Synthesize random initial state realization
    r0_err = randn(3,1)*simparams.sig_pos;
    v0_err = randn(3,1)*simparams.sig_vel;
    
    % r0_err = zeros(3,1);
    % v0_err = zeros(3,1);
    del_x_0 = [r0_err; v0_err]
    x0 = x_opt(1:6,1) + del_x_0;
    
    % Propagate initial error state thru trajectory
    event_times = [maneuver_times, tcm_target_time];
    event_dvs = [deltaVs_nom, zeros(3,1)];
    [x_tt, stm_tt, tt, x_currt] = mcProp(x0,event_times,event_dvs,simparams);
    
    
    
    % Use STM to correct del_r via tcm_dv executed at tcm_time
    tcm_idx(1) = find( tcm_time == tt );
    tcm_idx(2) = find( tcm_target_time == tt );

    % Determine tcm target dispersion
    del_x_n_num = x_tt(:,end) - tcm_target
    
    stmN0 = stm_tt(:,:,tcm_idx(2));
    J = [zeros(3,3), eye(3,3); -eye(3,3), zeros(3,3)];
    stmC0 = stm_tt(:,:,tcm_idx(1));
    stm0C = -J * stmC0' * J; % symplectic inverse
    stmNC = stmN0 * stm0C;
    
    T = [-inv( stmNC(1:3,4:6) ) * stmNC(1:3,1:3), -eye(3)];
    

    del_x_c = x_tt(:,tcm_idx(1)) - x_t(t==tcm_time,:)'
    tcm_dv = T * del_x_c


    tcm_norm_sav(end+1) = norm(tcm_dv);
    
    % test out stm linear propagation of del_x_0 to del_x_f
    del_x_n_lin = stmNC * stmC0 * del_x_0;
    linPropErrorPercent = (del_x_n_num - del_x_n_lin)./del_x_n_num * 100
    
    
    
    %% Propagate again with correction, see how close it comes
    event_times = [maneuver_times, tcm_time, tcm_target_time];
    event_dvs = [deltaVs_nom, tcm_dv, zeros(3,1)];
    [x_te, stm_te, te, x_curre] = mcProp(x0,event_times,event_dvs,simparams);
    
    del_xe = x_curre - tcm_target
    
    del_xe_f = stmNC * (del_x_0 + [zeros(3,1); tcm_dv])

    tcm_norm_sav(end) = tcm_norm_sav(end) + norm(del_xe(4:6));
    %% plot
    
    plot3(x_tt(1,:),x_tt(2,:),x_tt(3,:))
    
    
%     plot3(x_te(1,:),x_te(2,:),x_te(3,:))

end

threesigma_tcm = 3 * std(tcm_norm_sav)
tcm_min
tcm_min / threesigma_tcm

fctr = sqrt(pchisq(3,0.68));

tcm_min * fctr
%%%%% The difference is pretty large
