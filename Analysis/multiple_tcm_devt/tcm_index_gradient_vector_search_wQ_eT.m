function [TCMr_idx_best, TCMr_time_best, minDV] = tcm_index_gradient_vector_search_wQ_eT(x, traj, TCMr_idx_best, vel_disp_flag, deltaV, minDV, modRange, simparams)
%tcm_index_gradient_search Alters the time index of each individual TCM in
%each direction (earlier and later) to search for an improvement in the
%total delta V
%et = entire trajectory - unlike the original version that only optimized
%a single targeting portion at a time

% traj contains at least: x_t, t, t_s, stm_t

% % OPTIONAL INPUT: P_i
% if length(varargin) < 1
%     P_i = simparams.P_initial;
% else
%     P_i = varargin{1};
% end




% if length(TCMr_idx_best) == 5
%     output = 1;
%     counter = 1;
%     fileID = fopen('5_tcm_gradient_steps.txt','w');
% else
%     output = 0;
% end

maneuver_seg = simparams.maneuverSegments;
maneuver_idxs = nan*ones(1,length(maneuver_seg));
for i = 1:length(maneuver_seg)
    match = find(traj.t_s==maneuver_seg(i),1)-1;
    if ~isempty(match)
        maneuver_idxs(i) = match;
    end
end

% counter = 0;

TCMr_time_best = traj.t(TCMr_idx_best)';

[event_times, event_indicator] = define_events_v2(x(:), traj.t, TCMr_time_best, simparams);

P_constrained_events = find(event_indicator == 3 | event_indicator==0);
P_constrained_events = P_constrained_events(2:end);

% event_idx_logical = logical(sum(traj.t'==event_times', 1));    
% event_idxs = find(event_idx_logical);



improving = 1;

repeat_idx_counter = 0;

while improving
% for j = 1:nIter

    gradient_vector = zeros(1,length(TCMr_idx_best));
    TCMr_idx_test_save = zeros(length(TCMr_time_best)-1 + 1); % One more than the loop below to add the gradient method to the last element
    minDV_save = ones(1,length(TCMr_time_best)-1 + 1) * 1e15;


%     for i = 1:length(TCMr_time_best)
    for i = modRange(1):modRange(2)

        TCMr_time_test = TCMr_time_best;
        TCMr_idx_test = TCMr_idx_best;
        
    
        % Go through each TCM time and modify in each direction (forwards and
        % backwards in time) to see if the delta V total improves. Continue in
        % that direction until it stops improving, update the time in the
        % baseline, then go on to the next.
    
        % Modify the i-th element backwards
        improved = 0;

        % To prevent the gradient from searching backwards on top of a
        % previous TCM
        if i == 1
            compare = 0;
        else
            compare = TCMr_idx_test(i-1);
        end 


        if TCMr_idx_test(i)-1 > compare
            

            % is the TCM being moved on top of another maneuver?
            if isempty(find(TCMr_idx_test(i)-1 == maneuver_idxs))
                % Moving one index to the left for the TCM DV test
                    TCMr_idx_test(i) = TCMr_idx_test(i)-1;
                    TCMr_time_test(i) = traj.t(TCMr_idx_test(i));
        
                    % Testing with the updated TCM 
%                     [P_target, testDV] = calc_covariance_tcmdv_v2(x, t, t_s, stm_t, stm_t_i, TCMr_time_test, vel_disp_flag, deltaV, P_i, range, simparams); 
                    
                    if length(TCMr_idx_test) ~= length(unique(TCMr_idx_test))
                        testDV = minDV;
                    else


                        Q_k_km1 = calc_Q_events(traj, x, TCMr_time_test, simparams);
                        [~, testDV, ~, P_i_minus] = calc_covariance_wQ_tcmdv_v3(x(:), traj, TCMr_time_test, vel_disp_flag, deltaV, simparams.P_initial, Q_k_km1, simparams);

                        % Extract the P constrained nodes position dispersion RSS

                        pos_disp_fails = false;
                        for p = 1:length(P_constrained_events)
                            pos_disp_fails = logical(sqrt(trace(P_i_minus(1:3,1:3,P_constrained_events(p)))) > simparams.P_max_r | pos_disp_fails);
                        end



%                         [P_target, testDV] = calc_covariance_wQ_tcmdv(x, traj, TCMr_time_test, vel_disp_flag, deltaV, P_i, simparams); 
                    end
        
                    % Checking if it was cheaper
                    if testDV < minDV
%                         if i == length(TCMr_time_best) && sqrt(trace(P_target(1:3,1:3))) > simparams.P_max_r
                        if pos_disp_fails
                            % If the position dispersion constraint is
                            % violated as a result of the mod
                            pp=1;
                        else
                            TCMr_idx_test_save(i,:) = TCMr_idx_test;
                            minDV_save(i) = testDV;
            
                            gradient_vector(i) = -1;
            
                            improved = 1;
                        end
        
        
                    end  
            end
        end
        
         
        % if going backwards didn't result in an improvement, try going
        % forwards instead
    
        if improved == 0

            TCMr_idx_test = TCMr_idx_best;
            TCMr_time_test = TCMr_time_best;

            % Preventing the test from going forward on top of the next TCM 
            if i == length(gradient_vector)
                compare = 1e15;
            else
                compare = TCMr_idx_test(i+1);
            end

            if TCMr_idx_test(i)+1 < compare
                % is the TCM being moved on top of another maneuver or past the end of the trajectory?
                if isempty(find(TCMr_idx_test(i)+1 == maneuver_idxs)) && TCMr_idx_test(i) + 1 < length(traj.t)
                    % Modifying forward one index
                    TCMr_idx_test(i) = TCMr_idx_test(i)+1;
                    TCMr_time_test(i) = traj.t(TCMr_idx_test(i));
    
                    % Calculating the TCM DV of the test index
    %                 [P_target, testDV] = calc_covariance_tcmdv_v2(x, t, t_s, stm_t, stm_t_i, TCMr_time_test, vel_disp_flag, deltaV, P_i, range, simparams); 
                    
    
                    if length(TCMr_idx_test) ~= length(unique(TCMr_idx_test))
                        testDV = minDV;
                    else


                        Q_k_km1 = calc_Q_events(traj, x, TCMr_time_test, simparams);
                        [~, testDV, ~, P_i_minus] = calc_covariance_wQ_tcmdv_v3(x(:), traj, TCMr_time_test, vel_disp_flag, deltaV, simparams.P_initial, Q_k_km1, simparams);

                        % Extract the P constrained nodes position dispersion RSS

                        pos_disp_fails = false;
                        for p = 1:length(P_constrained_events)
                            pos_disp_fails = logical(sqrt(trace(P_i_minus(1:3,1:3,P_constrained_events(p)))) > simparams.P_max_r | pos_disp_fails);
                        end

%                         [P_target, testDV] = calc_covariance_wQ_tcmdv(x, traj, TCMr_time_test, vel_disp_flag, deltaV, P_i, simparams); 
                    end
        
                    % Comparing
                    if testDV < minDV
%                         if i == length(TCMr_time_best) && sqrt(trace(P_target(1:3,1:3))) > simparams.P_max_r
                        if pos_disp_fails
                            ppp=1;
                        else
                            TCMr_idx_test_save(i,:) = TCMr_idx_test;
                            minDV_save(i) = testDV;
        
                            gradient_vector(i) = 1;
                            minDV = testDV;
                        end
    
                    end

    
                end
            end
    
        end
    
    end

    

    if isempty(find(gradient_vector~=0)) || repeat_idx_counter == 3
        improving = 0;
    else
        % Taking a step in the gradient direction
        TCMr_idx_best = TCMr_idx_best + gradient_vector;
        TCMr_time_best = traj.t(TCMr_idx_best)';

        if length(TCMr_idx_best) ~= length(unique(TCMr_idx_best))
            % Do nothing
            repeat_idx_counter = repeat_idx_counter + 1;
        else

    
            % Storing a new min DV after the step
    %         [~, minDV] = calc_covariance_tcmdv_v2(x, t, t_s, stm_t, stm_t_i, TCMr_time_best, vel_disp_flag, deltaV, P_i, range, simparams); 
%             [~, minDV] = calc_covariance_wQ_tcmdv(x, traj, TCMr_time_best, vel_disp_flag, deltaV, P_i, simparams); 




            Q_k_km1 = calc_Q_events(traj, x, TCMr_time_best, simparams);
            [~, testDV, ~, P_i_minus] = calc_covariance_wQ_tcmdv_v3(x(:), traj, TCMr_time_test, vel_disp_flag, deltaV, simparams.P_initial, Q_k_km1, simparams);




            pos_disp_fails = false;
            for p = 1:length(P_constrained_events)
                pos_disp_fails = logical(sqrt(trace(P_i_minus(1:3,1:3,P_constrained_events(p)))) > simparams.P_max_r | pos_disp_fails);
            end


            if ~ pos_disp_fails
                TCMr_idx_test_save(end,:) = TCMr_idx_best;
                minDV_save(end) = minDV;
            end
        end

        % Check against the saved individual mods for the lowest
        [minDV, minIdx] = min(minDV_save);

        % Leaving this commented code here to output the search to a file.

%         if output
%             load('orbital_params.mat');
%             Rm = moon.a; 
%             n = sqrt((earth.mu + moon.mu)/Rm^3);
%             
%             ndTime2sec = 1/n;
%             ndTime2hrs = 1/n/3600;
%             ndTime2days = 1/n/3600/24;
%             ndDist2km = Rm;
%             ndVel2kms = Rm * n;
% 
%             outline = [counter, traj.t(TCMr_idx_test_save(minIdx,:))'*ndTime2hrs, minDV*ndVel2kms*1000];
% 
% 
%             fprintf(fileID,'%7.6f %7.6f %7.6f %7.6f %7.6f %7.6f\r\n',outline);
% 
% 
%             counter = counter + 1;
% 
% 
%         end



        TCMr_idx_best = TCMr_idx_test_save(minIdx,:);
        TCMr_time_best = traj.t(TCMr_idx_best)';


%         counter = counter + 1;

    end


end

% if output
%     fclose(fileID);
%     ppp=1;
% end


end