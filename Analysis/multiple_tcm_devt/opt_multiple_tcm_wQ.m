function [tcm_time_combined, tcm_idx_combined, min_tcm_dv_total, P_i_minus, P_i_plus, tcm_dv_each] = opt_multiple_tcm_wQ(x, traj, deltaVs_nom, simparams)
%opt_multiple_tcm Determines the optimal number of TCMs to perform along a
%trajectory, the times to perform them, and total 1 SIGMA TCM delta V.

% Note: working in time indices rather than actual times for the majority
% of what is below

%% Output variables
tcm_time_combined = [];
tcm_idx_combined = [];
min_tcm_dv_total = 0;
tcm_dv_each = [];
% event_logical = [];
% P_i_minus = [];
% P_i_plus = [];
% Removed outputs: tcm_time_cell,tcm_idx_cell,tcm_num_option_DVs


%% Setup

x = reshape(x,simparams.m,simparams.n);
P = simparams.P_initial;

vel_disp_flag = 0;

%% Loop through the position dispersion constraints stored in simparams.P_constrained_nodes 
% and calculate the optimal TCMs for each independent coast segment

for i = 1:length(simparams.P_constrained_nodes)
    
    if i == 1 && simparams.start_P_growth_node == 1

        start_idx = 1;

        target_node = simparams.P_constrained_nodes(i);
        target_time = sum(x(7,1:target_node - 1));
        target_idx = find(target_time == traj.t);


        assert(~isempty(target_idx));

        if length(target_idx) > 1
            target_idx = target_idx(1);
        end

        traj_eval.t = traj.t(start_idx:target_idx);
        traj_eval.t_s = traj.t_s(start_idx:target_idx);
        traj_eval.stm_t = traj.stm_t(:,:,start_idx:target_idx);

        % Attempt at fixing a bug: pass all of stm_t_i regardless of
        % segmentation
        traj_eval.stm_t_i = traj.stm_t_i(1:target_node - 1);
%         traj_eval.stm_t_i = traj.stm_t_i;
        
        
        
        traj_eval.x_t = traj.x_t(start_idx:target_idx,:);
%         traj_eval.x_i_f = traj.x_i_f(:,1:target_node-1);
%         traj_eval.stm_i = traj.stm_i(:,:,1:target_node-1);
        



        % Process noise portions
        % Process noise by time element accumulated from t0:
        traj_eval.Q_t = traj.Q_t(:,:,1:target_idx);
        % Process noise cells accumulated from the beginning of each seg:
        traj_eval.Q_t_i = traj.Q_t_i(1:target_node-1);

        

        
    else

        if i == 1
            start_time = sum(x(7,1:simparams.start_P_growth_node-1));
            start_idx = find(traj.t == start_time);
            start_node = simparams.start_P_growth_node;
        else
            start_idx = target_idx;
            start_node = target_node;
        end

        target_node = simparams.P_constrained_nodes(i);
        target_time = sum(x(7,1:target_node - 1));
        target_idx = find(target_time == traj.t, 1);

        assert(~isempty(target_idx));

        traj_eval.t = traj.t(start_idx:target_idx);
        traj_eval.t_s = traj.t_s(start_idx:target_idx);
        traj_eval.x_t = traj.x_t(start_idx:target_idx,:);

        % Want the STM from the beginning of the current correction portion (M)
        % to the end of the current correction portion (N): stmNtM
        % but currently have stmM0 and stmNt0
        % stmNtM = stmNt0 * stm0M

        stm0M = invert_stm(traj.stm_t(:,:,start_idx), simparams);
        traj_eval.stm_t = tmult(traj.stm_t(:,:,start_idx:target_idx), stm0M);


        traj_eval.stm_t_i = traj.stm_t_i(start_node:target_node-1);
%         traj_eval.stm_t_i = traj.stm_t_i;

        
%         traj_eval.stm_i = traj.stm_i(:,:,start_node:target_node-1);
%         traj_eval.x_i_f = traj.x_i_f(:,start_node:target_node-1);



        % Process noise by time, modified to accumulate from the beginning
        % of the new eval portion (where P was previously targeted):
        Qbarstart_tstart_t0 = traj.Q_t(:,:,start_node);
        Qbar_t_tstart = tmult(traj_eval.stm_t, tmult(Qbarstart_tstart_t0, traj_eval.stm_t,[0 1]));
        traj_eval.Q_t = traj.Q_t(:,:,start_idx:target_idx) - Qbar_t_tstart;

        % Process noise cells accumulated from the beginning of each seg,
        % corresponding to the eval segments:
        traj_eval.Q_t_i = traj.Q_t_i(start_node:target_node-1);

    end

    % Get all deltaVs that occur at start_idx and less than target_idx
    maneuver_include = simparams.maneuverSegments>=traj.t_s(start_idx) & simparams.maneuverSegments < target_node;
    deltaVs_nom_eval = deltaVs_nom(:,maneuver_include);
    
    % Only include the final velocity dispersion correction if it is the
    % final target node
    if i == length(simparams.P_constrained_nodes)
        vel_disp_flag = 1;
        deltaVs_nom_eval = [deltaVs_nom_eval, deltaVs_nom(:,end)];
    end
    
    %% First, calculate the location of the lowest Delta V option that meets the
    % target position dispersion constraint
    
    % Structure for storing the delta Vs: tcm_num_option_DVs - each
    % index corresponds to the total DV for the number of corresponding TCMrs
%     [tcm_time,tcm_idx,tcm_num_option_DVs] = min_dv_tcm_meets_dispersion_constraint(x(:), t_eval, t_s_eval, stm_t_eval, stm_t_i_eval, vel_disp_flag, deltaVs_nom_eval, P, simparams);
    range = [start_idx, target_idx];
%     [tcm_time,tcm_idx,tcm_num_option_DVs] = min_dv_tcm_meets_dispersion_constraint_v2(x(:), t, t_s, stm_t, stm_t_i, vel_disp_flag, deltaVs_nom_eval, P, range, simparams);



    [tcm_time,tcm_idx,tcm_num_option_DVs] = min_dv_tcm_meets_dispersion_constraint_wQ_v2(x(:), traj_eval, vel_disp_flag, deltaVs_nom_eval, P, simparams);


    
    %% Find the optimal number and time of TCMs
    [tcm_time,tcm_idx,min_tcm_dv] = optimize_tcm_guess_wQ(x(:), traj_eval, tcm_time, tcm_idx, tcm_num_option_DVs, vel_disp_flag, deltaVs_nom_eval, P, simparams);
%     [tcm_time,tcm_idx,min_tcm_dv] = optimize_tcm_guess(x(:), t, t_s, stm_t, stm_t_i, tcm_time, tcm_idx, tcm_num_option_DVs, vel_disp_flag, deltaVs_nom_eval, P, range, simparams);
    P = calc_covariance_wQ_tcmdv(x(:), traj_eval, tcm_time, vel_disp_flag, deltaVs_nom_eval, P, simparams);
   



    % Calculate the covariance at the end of the coast portion to use at
    % the beginning of the next TCM coast portion
    Pprev=P;


%     [P, tcm_dv_total, tcm_dv, P_i_minus_portion, P_i_plus_portion] = calc_covariance_tcmdv_v2(x(:), t_eval, t_s_eval, stm_t_eval, stm_t_i_eval, tcm_time, vel_disp_flag, deltaVs_nom_eval, P, range, simparams);
%     [P, tcm_dv_total, tcm_dv, P_i_minus_portion, P_i_plus_portion] = calc_covariance_tcmdv_v2(x(:), t, t_s, stm_t, stm_t_i, tcm_time, vel_disp_flag, deltaVs_nom_eval, P, range, simparams);

    if isnan(P(1))
        ppp=1;
    end


%     if i == 1
%         P_i_minus = P_i_minus_portion;
%         P_i_plus = P_i_plus_portion;
%     else
%         P_i_minus(:,:,end+1:end+size(P_i_minus_portion,3)) = P_i_minus_portion;
%         P_i_plus(:,:,end+1:end+size(P_i_plus_portion,3)) = P_i_plus_portion;
%     end

    %% Store in output data structures
    tcm_time_combined = [tcm_time_combined, tcm_time];
    if range(1) > 1
        tcm_idx_combined = [tcm_idx_combined, tcm_idx + start_idx - 1];
    else
        tcm_idx_combined = [tcm_idx_combined, tcm_idx];
    end
    min_tcm_dv_total = min_tcm_dv_total + min_tcm_dv;

%     tcm_dv_each = [tcm_dv_each, tcm_dv(1:end-1+vel_disp_flag)];
%     tcm_dv_each = [tcm_dv_each, tcm_dv];

end

% Incorporate the cost of the second portion when optimizing the first
% portion
% modRange = [1, length(tcm_time_combined) - length(tcm_time) + 1];
% [TCMr_idx_best, TCMr_time_best, minDV] = tcm_index_gradient_vector_search_wQ_eT(x, traj, tcm_idx_combined, vel_disp_flag, deltaVs_nom, min_tcm_dv_total, modRange, simparams);
% tcm_idx_combined = TCMr_idx_best;
% tcm_time_combined = TCMr_time_best;

% Recalculate with more accurate function (shorter STM segment portions vs
% reaching back to the beginning of the trajectory like
% calc_covariance_tcmdv does above)

Q_k_km1 = calc_Q_events(traj, x, tcm_time_combined, simparams);

[P, min_tcm_dv_total, tcm_dv_each, P_i_minus, P_i_plus] = calc_covariance_wQ_tcmdv_v3(x(:), traj, tcm_time_combined, vel_disp_flag, deltaVs_nom, simparams.P_initial, Q_k_km1, simparams);




end