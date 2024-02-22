function [J, J_gradient] = obj_min_tcm(x, simparams)
%obj_min_tcm Calculates the objective function value
%   The objective function value is the sum of the nominal impulsive delta
%   V plus two estimated 3 sigma correction delta V values.



m = simparams.m; % the number of elements per segment
n = simparams.n; % the number of segments
nsv = simparams.nsv; % number of state variables

% Reshaping x so each column is a segment initial state and duration
x = reshape(x,m,n);

% The three parameters that form the objective function
% deltaV, dvR3sigma, and dvV3sigma

% If want to output gradients, this helps:
if nargout > 1
    outputGradients = 1;
else
    outputGradients = 0;
end


% NO GRADIENT OUTPUT
if ~outputGradients
 
    if simparams.perform_correction
        % Propagate traj state, STM, and Q history
        traj = createStateStmQHistory(x, simparams);
        % Using saved dynamics, calculate total impulsive delta V
        [deltaV, deltaVs_nom] = calcDeltaV(x, traj.x_i_f, traj.stm_i, simparams);
        % Function to find the optimal number of TCMs and times to perform them along the nominal trajectory
%         [~, ~, min_tcm_dv] = opt_multiple_tcm_wQ(x, traj, deltaVs_nom, simparams);
        % TCM at nodes method
        tcm_time = zeros(1, length(simparams.tcm_nodes));
%         tcm_idx = zeros(1, length(simparams.tcm_nodes));

        for i = 1:length(simparams.tcm_nodes)                
            tcm_time(i) = sum(x(m,1:simparams.tcm_nodes(i)-1));
        end
%         for i = 1:length(tcm_time)
%             tcm_idx(i) = find(traj.t == tcm_time(i))';
%         end

        % Calculate the process noise covariance between updates
        [Q_k_km1] = calc_Q_events(traj, x, tcm_time, simparams);
        % Calculate the TCM dv and dispersion covariance prior to each covariance modifying event (k)            
        [~, min_tcm_dv] = calc_covariance_wQ_tcmdv_v3(x(:), traj, tcm_time, 1, deltaVs_nom, simparams.P_initial, Q_k_km1, simparams);

        tcm_3sigma = simparams.tcm_rss_factor * min_tcm_dv;
    else
        traj = createStateStmHistory(x, simparams);
        deltaV = calcDeltaV(x, traj.x_i_f, traj.stm_i, simparams);
        tcm_3sigma = 0;
    end

    
% YES GRADIENT OUTPUT
else
    % Propagate entire trajectory and necessary things for gradients
    if simparams.perform_correction
        if sum(sum(simparams.Qt)) == 0 % If there is no process noise, don't prop Qbar and dQbar
            traj  = createStateStmSttHistory(x, simparams);
        else
            traj = createStateStmSttQdQHistory(x, simparams);
        end        
    else
        traj = createStateStmHistory(x, simparams);    
    end

    if existsAndTrue('rdvz_flag', simparams)
        x = reshape(x, simparams.m, simparams.n);
        totalTime = sum(x(m,:));
        simparams.x_target = stateProp(simparams.x0_target, totalTime, simparams);
    end

    % Using saved dynamics, calculate total nominal impulsive delta V and
    % the delta V gradient
    [deltaV, deltaVs_nom, deltaV_gradient] = calcDeltaV(x, traj.x_i_f, traj.stm_i, simparams);
%     tst = calc_deltaV_gradient(x, x_i_f, stm_i, simparams);

    % Using saved dynamics, calculate the 3 sigma TCM pair

    if simparams.perform_correction
         
        
        
          
    
        % Calculate the TCM gradient
        
        if sum(sum(simparams.Qt)) == 0
            % Function to optimize the number and location of TCMs along the trajectory
%             [tcm_time, tcm_idx, min_tcm_dv, P_i_minus, P_i_plus] = opt_multiple_tcm(x, deltaVs_nom, traj.t, traj.t_s, traj.stm_t, traj.stm_t_i, simparams);

            tcm_time = zeros(1, length(simparams.tcm_nodes));
            tcm_idx = zeros(1, length(simparams.tcm_nodes));

            for i = 1:length(simparams.tcm_nodes)                
                tcm_time(i) = sum(x(m,1:simparams.tcm_nodes(i)-1));
            end
            for i = 1:length(tcm_time)
                tcm_idx(i) = find(traj.t == tcm_time(i))';
            end

            tcm_gradient = calc_multiple_tcm_gradient(x, traj.x_i_f, traj.stm_i, traj.stm_t, traj.stm_t_i, traj.stt_t_i, traj.t, traj.t_s, tcm_time, tcm_idx, P_i_minus, deltaVs_nom, simparams);
        else
            % Optimize the number and location of TCMs with process noise
%             [tcm_time, tcm_idx, min_tcm_dv, P_i_minus] = opt_multiple_tcm_wQ(x, traj, deltaVs_nom, simparams);
            % TCM at nodes method
            tcm_time = zeros(1, length(simparams.tcm_nodes));
            tcm_idx = zeros(1, length(simparams.tcm_nodes));

            for i = 1:length(simparams.tcm_nodes)                
                tcm_time(i) = sum(x(m,1:simparams.tcm_nodes(i)-1));
            end
            for i = 1:length(tcm_time)
                tcm_idx(i) = find(traj.t == tcm_time(i))';
            end


            % Calculate the process noise and process noise sensitivity tensors
            [Q_k_km1, dQ_k_km1_dxi, dQ_k_km1_ddti] = calc_Q_events(traj, x, tcm_time, simparams);
            % Calculate the TCM dv and dispersion covariance prior to each covariance modifying event (k)            
            [~, min_tcm_dv, ~, P_i_minus] = calc_covariance_wQ_tcmdv_v3(x(:), traj, tcm_time, 1, deltaVs_nom, simparams.P_initial, Q_k_km1, simparams);


            
            tcm_gradient = calc_multiple_tcm_gradient_wQ(x, traj, tcm_time, tcm_idx, P_i_minus, dQ_k_km1_dxi, dQ_k_km1_ddti, deltaVs_nom, simparams);
        end

        % Multiply the TCM RSS by the factor (3 sigma mostly)
        tcm_3sigma = simparams.tcm_rss_factor * min_tcm_dv;
        
        tcm_gradient = simparams.tcm_rss_factor * tcm_gradient;
    else
        tcm_3sigma = 0;
        tcm_gradient = 0 * x(:);
    end



    


end



%% 

% Objective function sum
J = deltaV + tcm_3sigma;
% J = deltaV;

% J = tcm_3sigma; %%%% DEBUGGING

% Gradient sum
if outputGradients
    J_gradient = deltaV_gradient' + tcm_gradient;
%     J_gradient = deltaV_gradient';
%     J_gradient = tcm_gradient; %%%% DEBUGGING
end




end