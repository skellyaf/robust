function [J, J_gradient] = obj_min_tcm(x, simparams)
%obj_min_tcm Calculates the objective function value
%   The objective function value is the sum of the nominal impulsive delta
%   V plus two estimated 3 sigma correction delta V values.



% m = simparams.m; % the number of elements per segment
% n = simparams.n; % the number of segments

% Reshaping x so each column is a segment initial state and duration
% x = reshape(x,m,n);

% The three parameters that form the objective function
% deltaV, dvR3sigma, and dvV3sigma

% If want to output gradients, this helps:
if nargout > 1
    outputGradients = 1;
else
    outputGradients = 0;
end



if ~outputGradients
    % Propagate entire trajectory state and STM history; save dynamics at each time step
    [stm_i, x_i_f, x_t, stm_t, t, t_s] = createStateStmHistory(x, simparams);

    % Using saved dynamics, calculate total impulsive delta V
    [deltaV, deltaVs_nom] = calcDeltaV(x, x_i_f, simparams);

    % Using saved dynamics, calculate the 3 sigma TCM pair

    if simparams.perform_correction
        % Function to find the minimum pair along the trajectory
%         [tcm_3sigma,tcm_time, dvR3sigma_tr, dvV3sigma_tr] = tcmPair_rv(x, t, stm_t, deltaVs_nom, simparams);
        % Function to find the optimal number of TCMs and times to perform
        % them along the nominal trajectory
        [~,~,min_tcm_dv] = opt_multiple_tcm(x, t, t_s, stm_t, simparams);
        tcm_3sigma = 3*min_tcm_dv;
    else
        tcm_3sigma = 0;
    end

else
    % Propagate entire trajectory state, STM and STT history; save dynamics at each time step
    % Need STTs for the TCM gradients    
    if simparams.perform_correction
        [stm_i, stt_i, x_i_f, x_t, stm_t, stt_t_i, t, t_s, stm_t_i]  = createStateStmSttHistory(x, simparams);
    else
        [stm_i, x_i_f, x_t, stm_t, t, t_s] = createStateStmHistory(x, simparams);

    end

    % Using saved dynamics, calculate total nominal impulsive delta V
    [deltaV, deltaVs_nom] = calcDeltaV(x, x_i_f, simparams);

    % Using saved dynamics, calculate the 3 sigma TCM pair

    if simparams.perform_correction
        % Function to optimize the number and location of TCMs along the trajectory 
        [tcm_time, tcm_idx, min_tcm_dv, P_i_minus, P_i_plus] = opt_multiple_tcm(x, deltaVs_nom, t, t_s, stm_t, stm_t_i, simparams);
        
        % Multiply the TCM RSS by 3
        tcm_3sigma = 3*min_tcm_dv;

        % Calculate the process noise and process noise sensitivity tensors
        [Q_k_km1, dQ_k_km1_dxi, dQ_k_km1_ddti] = calc_Q_events(traj, x, tcm_time, simparams);
    
        % Calculate the TCM gradient
        tcm_gradient = calc_multiple_tcm_gradient_wQ(x, x_t, x_i_f, stm_i, stt_i, stm_t, stm_t_i, stt_t_i, t, t_s, tcm_time, tcm_idx, P_i_minus, dQ_k_km1_dxi, dQ_k_km1_ddti, deltaVs_nom, simparams);
        tcm_gradient = 3*tcm_gradient;
    else
        tcm_3sigma = 0;
        tcm_gradient = 0 * x(:);
    end


    % Calculate objective function gradient
    % Delta V gradient
    deltaV_gradient = calc_deltaV_gradient(x, x_i_f, stm_i, simparams);


end



%% 

% Objective function sum
J = deltaV + tcm_3sigma;

% J = tcm_3sigma; %%%% DEBUGGING

% Gradient sum
if outputGradients
    J_gradient = deltaV_gradient' + tcm_gradient;
%     J_gradient = tcm_gradient; %%%% DEBUGGING
end




end