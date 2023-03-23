function [tcm_time,tcm_idx,min_tcm_dv,tcm_time_cell,tcm_idx_cell,tcm_num_option_DVs] = opt_multiple_tcm(x, t, stm_t, simparams)
%opt_multiple_tcm Determines the optimal number of TCMs to perform along a
%trajectory, the times to perform them, and total 1 SIGMA TCM delta V.

% Description

% Note: working in time indices rather than actual times for the majority
% of what is below


%% Setup

x = reshape(x,simparams.m,simparams.n);

% Determine target info
% STM from beginning of trajectory to the state being targeted

if simparams.target_final_maneuver
    final_seg = simparams.maneuverSegments(end);
    target_time = sum(x(7,1:final_seg - 1));
    target_idx = find(target_time == t) ;
    stmN0 = stm_t(:,:,target_idx);
    % Abbreviate t and stm_t to go only until the target (not the end of
    % the trajectory)
    t = t(1:target_idx);
    stm_t = stm_t(:,:,1:target_idx);
else
    stmN0 = stm_t(:,:,end);    
end



%% First, calculate the location of the lowest Delta V option that meets the
% target position dispersion constraint

% Structure for storing the delta Vs: tcm_num_option_DVs - each
% index corresponds to the total DV for the number of corresponding TCMrs

[tcm_time,tcm_idx,tcm_num_option_DVs] = min_dv_tcm_meets_dispersion_constraint(x, t, stm_t, simparams);

%% Structure to store the TCM info
tcm_time_cell{1} = tcm_time; % Cell structure
tcm_idx_cell{1} = tcm_idx;


%% Incrementally loop through adding a maneuver, optimizing, and storing the Delta V
% until the Delta V stops decreasing

num_tcmr = length(tcm_idx);
improving = 1;
while improving

    tcm_time = tcm_time_cell{end};
    tcm_idx = tcm_idx_cell{end};

    % Add a maneuver at the lowest delta V time option
    
    % Create a logical array of indices to test
%     test_logical = logical(ones(1,length(t)));
    test_logical = true(1,length(t));

    if length(tcm_idx) == 1
        test_logical(tcm_idx:end) = false;
    else
        test_logical(tcm_idx) = false;
        test_logical(1:tcm_idx(1)) = false;
        test_logical(tcm_idx(end)+1:end) = false;
    end

    plusOneTCMr_totalDV_t = ones(1,length(t)) * nan;
%     plusOneTCMr_totalDV_t = ones(1,length(t)) * 1e8;

    % Test each index in the logical array

    if sum(test_logical)==0
        [~, min_tcm_dv] = calc_covariance_tcmdv(x, t, stm_t, tcm_time, simparams);

        improving = 0;
    else
        for i = find(test_logical)
            tcm_new_i = t(i);
            tcm_time_test = sort([tcm_time, tcm_new_i]);
            [~, plusOneTCMr_totalDV_t(i)] = calc_covariance_tcmdv(x, t, stm_t, tcm_time_test, simparams);         
        end
       

        % Find the minimum
        [minDV, minIdx] = min(plusOneTCMr_totalDV_t);
        % Add the corresponding minimum time to a test vector
        tcm_time_test = sort([tcm_time, t(minIdx)]);
        tcm_idx_test = sort([tcm_idx, minIdx]);
    
        % Perform gradient descent on the indices

        if ~isempty(find(tcm_idx_test(1) == tcm_idx_test(2:end)))
            ppp=1;
        end

        [tcm_idx_test, tcm_time_test, minDV] = tcm_index_gradient_vector_search(x, t, stm_t, tcm_idx_test, 1, simparams);
        
        if length(tcm_idx_test) > 2
            % Perform stochastic gradient descent to get double indices
            [tcm_idx_test] = random_unit_search(x, t, stm_t, tcm_idx_test, 1, simparams);
        
            % Perform gradient descent one more time on the indices
            [tcm_idx_test, tcm_time_test, minDV] = tcm_index_gradient_vector_search(x, t, stm_t, tcm_idx_test, 1, simparams);
        end
    
        % Save it as the minimum:
        % Save the TCM times in a structure
        tcm_time_cell{end+1} = tcm_time_test;
        % Save the TCM indices in a structure
        tcm_idx_cell{end+1} = tcm_idx_test;
        % Save the total TCM delta V in a structure
        tcm_num_option_DVs(end+1) = minDV;
    
    
        % Check if the Delta V improved, end the while loop if not
        if minDV >= tcm_num_option_DVs(end-1)
            tcm_time = tcm_time_cell{end-1};
            tcm_idx = tcm_idx_cell{end-1};
            min_tcm_dv = tcm_num_option_DVs(end-1); %%%% 1 SIGMA ONLY
            improving = 0;
        end 

        if length(tcm_idx_test) == 6 %%%%% MAXIMUM NUMBER OF TCMS ALLOWED
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