function [tcm_time_combined, tcm_idx_combined, min_tcm_dv_total, P_i_minus, P_i_plus, tcm_dv_each] = opt_multiple_tcm(x, t, t_s, stm_t, simparams)
%opt_multiple_tcm Determines the optimal number of TCMs to perform along a
%trajectory, the times to perform them, and total 1 SIGMA TCM delta V.
%% MODIFICATION 23 MAR 23: IMPLEMENTING A MID-TRAJECTORY POSITION DISPERSION CONSTRAINT




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
    
    if i == 1
        start_idx = 1;
        target_node = simparams.P_constrained_nodes(i);
        target_time = sum(x(7,1:target_node - 1));
        target_idx = find(target_time == t);

        assert(length(target_idx)>0);
        if length(target_idx) > 1
            target_idx = target_idx(1);
        end

        t_eval = t(start_idx:target_idx);
        t_s_eval = t_s(start_idx:target_idx);
        stm_t_eval = stm_t(:,:,start_idx:target_idx);
        
    else

        start_idx = target_idx;

        target_node = simparams.P_constrained_nodes(i);
        target_time = sum(x(7,1:target_node - 1));
        target_idx = find(target_time == t);

        assert(length(target_idx)>0);

        if length(target_idx) > 1
            target_idx = target_idx(1);
        end

        t_eval = t(start_idx:target_idx);
        t_s_eval = t_s(start_idx:target_idx);

        % Want the STM from the beginning of the current correction portion (M)
        % to the end of the current correction portion (N): stmNtM
        % but currently have stmM0 and stmNt0
        % stmNtM = stmNt0 * stm0M
        testStm = stm_t(:,:,start_idx);
        if size(testStm,3) ~= 1
            ppp=1; % debug
        end
        stm0M = invert_stm(stm_t(:,:,start_idx), simparams);
        stm_t_eval = tmult(stm_t(:,:,start_idx:target_idx), stm0M);

    end
    
    % Only include the final velocity dispersion correction if it is the
    % final target node
    if i == length(simparams.P_constrained_nodes)
        vel_disp_flag = 1;
    end
    
    %% First, calculate the location of the lowest Delta V option that meets the
    % target position dispersion constraint
    
    % Structure for storing the delta Vs: tcm_num_option_DVs - each
    % index corresponds to the total DV for the number of corresponding TCMrs
    [tcm_time,tcm_idx,tcm_num_option_DVs] = min_dv_tcm_meets_dispersion_constraint(x(:), t_eval, t_s_eval, stm_t_eval, vel_disp_flag, P, simparams);
    
    %% Find the optimal number and time of TCMs
    [tcm_time,tcm_idx,min_tcm_dv] = optimize_tcm_guess(x(:), t_eval, t_s_eval, stm_t_eval, tcm_time, tcm_idx, tcm_num_option_DVs, vel_disp_flag, P, simparams);

    % Calculate the covariance at the end of the coast portion to use at
    % the beginning of the next TCM coast portion


%     if i == 2
%         tcm_time(1) = sum(x(7,1:11));
%     end

    [P, tcm_dv_total, tcm_dv, P_i_minus_portion, P_i_plus_portion] = calc_covariance_tcmdv(x(:), t_eval, t_s_eval, stm_t_eval, tcm_time, vel_disp_flag, P, simparams);


    if i == 1
        P_i_minus = P_i_minus_portion;
        P_i_plus = P_i_plus_portion;
    else
        P_i_minus(:,:,end+1:end+size(P_i_minus_portion,3)) = P_i_minus_portion;
        P_i_plus(:,:,end+1:end+size(P_i_plus_portion,3)) = P_i_plus_portion;
    end

    %% Store in output data structures
    tcm_time_combined = [tcm_time_combined, tcm_time];
    tcm_idx_combined = [tcm_idx_combined, tcm_idx + start_idx - 1];
    min_tcm_dv_total = min_tcm_dv_total + min_tcm_dv;

%     tcm_dv_each = [tcm_dv_each, tcm_dv(1:end-1+vel_disp_flag)];
    tcm_dv_each = [tcm_dv_each, tcm_dv];

end




end