function [tcm_time,tcm_idx,min_tcm_dv,tcm_time_cell,tcm_idx_cell,tcm_num_option_DVs] = optimize_tcm_guess_wQ(x, traj, tcm_time, tcm_idx, tcm_num_option_DVs, vel_disp_flag, deltaV, P_i, simparams)
%optimize_tcm_guess takes the single final TCM to meet the position
%dispersion constraint and incrementally adds a TCM and optimizes the TCMs
%until the lowest TCM magnitude is found (stops searching when adding a TCM
%increases the total cost)

% OPTIONAL INPUT: P_i
% if length(varargin) < 1
%     P_i = simparams.P_initial;
% else
%     P_i = varargin{1};
% end


%% Structure to store the TCM info
tcm_time_cell{1} = tcm_time; % Cell structure
tcm_idx_cell{1} = tcm_idx;


%% Incrementally loop through adding a maneuver, optimizing, and storing the Delta V
% until the Delta V stops decreasing

% num_tcmr = length(tcm_idx);
improving = 1;
while improving

    
    tcm_time = tcm_time_cell{end};
    tcm_idx = tcm_idx_cell{end};

    % Add a maneuver at the lowest delta V time option
    
    % Create a logical array of indices to test
    test_step_size = 10; %%%%%%%%% INSTEAD OF SEARCHING EVERY INDEX ALONG THE LINE, ONLY SEARCH EVERY X INDICES, THEN GRADIENT DESCEND TO THE MIN IN THE NEXT STEP 
    
    test_logical = false(1,length(traj.t));
    test_logical(1:test_step_size:end-1) = true;
    test_logical(1:25) = true;

%     test_logical(range(1):test_step_size:range(2)-1) = true;
%     test_logical(range(1):range(1)+25) = true;

    % Don't test where nominal maneuvers are if we're already correcting
    % the nominal maneuvers there
    if simparams.correct_nominal_dvs
        for i = 1:length(simparams.maneuverSegments)
            maneuverSegIdx = traj.t_s==simparams.maneuverSegments(i);
            maneuverIdx = find(maneuverSegIdx,1)-1; % the time of the maneuver occurs at the intersection of 1 and 2, and the previous segment gets the final time index
            test_logical(maneuverIdx) = false;
        end
    end




    if length(tcm_idx) == 1
        test_logical(tcm_idx:end) = false;
    else
        test_logical(tcm_idx) = false; % Don't test on top of other TCMs
        test_logical(1:tcm_idx(1)) = false; % Don't test before the first TCM
        test_logical(tcm_idx(end)+1:end) = false; % Don't test after the last TCM
    end


    plusOneTCMr_totalDV_t = ones(1,length(traj.t)) * nan;

    % Test each index in the logical array

    if sum(test_logical)==0
        [~, min_tcm_dv] = calc_covariance_wQ_tcmdv(x, traj, tcm_time, vel_disp_flag, deltaV, P_i, simparams);

        improving = 0;
    else
        
        for i = find(test_logical)
            tcm_new_i = traj.t(i);
            tcm_time_test = sort([tcm_time, tcm_new_i]);
%             [~, plusOneTCMr_totalDV_t(i)] = calc_covariance_tcmdv_v2(x, t, t_s, stm_t, stm_t_i, tcm_time_test, vel_disp_flag, deltaV, P_i, range, simparams);         
            [~, plusOneTCMr_totalDV_t(i)] = calc_covariance_wQ_tcmdv(x, traj, tcm_time_test, vel_disp_flag, deltaV, P_i, simparams);         
        end
        

        
       

        % Find the minimum...plot command: plot(traj.t(~isnan(plusOneTCMr_totalDV_t)), plusOneTCMr_totalDV_t(~isnan(plusOneTCMr_totalDV_t)))
        [minDV, minIdx] = min(plusOneTCMr_totalDV_t);
        % Add the corresponding minimum time to a test vector
        tcm_time_test = sort([tcm_time, traj.t(minIdx)]);
        tcm_idx_test = sort([tcm_idx, minIdx]);
    
        % Perform gradient descent on the indices

        [tcm_idx_test, tcm_time_test, minDV] = tcm_index_gradient_vector_search_wQ(x, traj, tcm_idx_test, vel_disp_flag, deltaV, P_i, minDV, simparams);

        if length(tcm_idx_test) > 2
            % Perform stochastic gradient descent to get double indices
            [tcm_idx_test] = random_unit_search_wQ(x, traj, tcm_idx_test, vel_disp_flag, deltaV, P_i, minDV, simparams);

            % Perform gradient descent one more time on the indices
            [tcm_idx_test, tcm_time_test, minDV] = tcm_index_gradient_vector_search_wQ(x, traj, tcm_idx_test, vel_disp_flag, deltaV, P_i, minDV, simparams);
        end
    


        % Save it as the minimum:
        % Save the TCM times in a structure
        tcm_time_cell{end+1} = tcm_time_test;
        % Save the TCM indices in a structure
        tcm_idx_cell{end+1} = tcm_idx_test;
        % Save the total TCM delta V in a structure
        tcm_num_option_DVs(end+1) = minDV;


    
    
        % Check if the Delta V (plus an improvement threshold) improved, end the while loop if not
        if minDV + simparams.add_tcm_improvement_threshold >= tcm_num_option_DVs(end-1)
            tcm_time = tcm_time_cell{end-1};
            tcm_idx = tcm_idx_cell{end-1};
            min_tcm_dv = tcm_num_option_DVs(end-1); %%%% 1 SIGMA ONLY
            improving = 0;
        

        elseif length(tcm_idx_test) == simparams.max_num_TCMs %%%%% MAXIMUM NUMBER OF TCMS ALLOWED
            tcm_time = tcm_time_cell{end};
            tcm_idx = tcm_idx_cell{end};
            min_tcm_dv = tcm_num_option_DVs(end); %%%% 1 SIGMA ONLY
            improving = 0;
        end

    end



end

% Return the lowest delta V tcm_time and tcm_idx arrays, the corresponding
% delta V

% Also return the structures containing the other number of TCM results


end