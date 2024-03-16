function [tcm_time, tcm_idx, min_tcm_dv_total, P_i_minus, P_i_plus, tcm_dv_each] = opt_multiple_tcm_wQ_multiPart(x, traj, deltaVs_nom, simparams)
%opt_multiple_tcm Determines the optimal number of TCMs to perform along a
%trajectory, the times to perform them, and total 1 SIGMA TCM delta V.

% Note: working in time indices rather than actual times for the majority
% of what is below

%% Setup

x = reshape(x,simparams.m,simparams.n);
P = simparams.P_initial;

vel_disp_flag = 0;

%% Loop through the position dispersion constraints stored in simparams.P_constrained_nodes 
% and calculate the optimal TCMs for each independent coast segment


% Step 1: add a single TCM for each target that is the cheapest way to meet the
% position dispersion constraint:
tcm_idx = [];
for i = 1:length(simparams.P_constrained_nodes)
    if i == length(simparams.P_constrained_nodes)
        vel_disp_flag = 1;
    end

    target_node = simparams.P_constrained_nodes(i);
    target_time = sum(x(7,1:target_node - 1));
    target_idx = find(target_time == traj.t);


    [tcm_time,tcm_idx,tcm_num_option_DVs, P_tgt_i] = min_dv_tcm_meets_dispersion_constraint_wQ_v3(x(:), traj, vel_disp_flag, deltaVs_nom, P, tcm_idx, target_idx, target_time, simparams);
    P = P_tgt_i;
end


% Step 2: Optimize the guess by adding TCMs at the cheapest spot and
% subsequently adjusting
tic
[tcm_time,tcm_idx,min_tcm_dv] = optimize_tcm_guess_wQ_entireTraj(x(:), traj, tcm_time, tcm_idx, tcm_num_option_DVs, vel_disp_flag, deltaVs_nom, simparams.P_initial, simparams);
toc

Q_k_km1 = calc_Q_events(traj, x, tcm_time, simparams);

[P, min_tcm_dv_total, tcm_dv_each, P_i_minus, P_i_plus] = calc_covariance_wQ_tcmdv_v3(x(:), traj, tcm_time, vel_disp_flag, deltaVs_nom, simparams.P_initial, Q_k_km1, simparams);


% [P, min_tcm_dv_total, tcm_dv_each, P_i_minus, P_i_plus]
%%%%%%%%%%%%%%%% COME BACK HERE, JUST FINISHED PLACING THE FINAL TCM OF
%%%%%%%%%%%%%%%% EACH SEGMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end