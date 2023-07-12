function [TCMr_idx_best, TCMr_time_best, minDV] = tcm_index_gradient_vector_search(x, t, t_s, stm_t, TCMr_idx_best, vel_disp_flag, deltaV, P_i, range, minDV, simparams)
%tcm_index_gradient_search Alters the time index of each individual TCM in
%each direction (earlier and later) to search for an improvement in the
%total delta V

% % OPTIONAL INPUT: P_i
% if length(varargin) < 1
%     P_i = simparams.P_initial;
% else
%     P_i = varargin{1};
% end




% if length(TCMr_idx_best) == 4
%     output = 1;
%     counter = 1;
%     fileID = fopen('4_tcm_gradient_steps.txt','w');
% else
%     output = 0;
% end

counter = 0;

TCMr_time_best = t(TCMr_idx_best)';

improving = 1;

% [~, minDV] = calc_covariance_tcmdv_v2(x, t, t_s, stm_t, stm_t_i, TCMr_time_best, vel_disp_flag, deltaV, P_i, range, simparams); 

while improving
% for j = 1:nIter

    gradient_vector = zeros(1,length(TCMr_idx_best));
    TCMr_idx_test_save = zeros(length(TCMr_time_best)-1 + 1); % One more than the loop below to add the gradient method to the last element
    minDV_save = ones(1,length(TCMr_time_best)-1 + 1) * 1e8;


    for i = 1:length(TCMr_time_best)

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
            % Moving one index to the left for the TCM DV test
            TCMr_idx_test(i) = TCMr_idx_test(i)-1;
            TCMr_time_test(i) = t(TCMr_idx_test(i));

            % Testing with the updated TCM 
%             [P_target, testDV] = calc_covariance_tcmdv_v2(x, t, t_s, stm_t, stm_t_i, TCMr_time_test, vel_disp_flag, deltaV, P_i, range, simparams); 
            [P_target, testDV] = calc_covariance_tcmdv(x, t, t_s, stm_t, TCMr_time_test, vel_disp_flag, deltaV, P_i, range, simparams); 

            % Checking if it was cheaper
            if testDV < minDV
                if i == length(TCMr_time_best) && sqrt(trace(P_target(1:3,1:3))) > simparams.P_max_r
                    % If it is the final element and the covariance
                    % constraint is violated

                else
                    TCMr_idx_test_save(i,:) = TCMr_idx_test;
                    minDV_save(i) = testDV;
    
                    gradient_vector(i) = -1;
    
                    improved = 1;
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
                compare = 1e8;
            else
                compare = TCMr_idx_test(i+1);
            end

            if TCMr_idx_test(i)+1 < compare
                % Modifying forward one index
                TCMr_idx_test(i) = TCMr_idx_test(i)+1;
                TCMr_time_test(i) = t(TCMr_idx_test(i));

                % Calculating the TCM DV of the test index
%                 [P_target, testDV] = calc_covariance_tcmdv_v2(x, t, t_s, stm_t, stm_t_i, TCMr_time_test, vel_disp_flag, deltaV, P_i, range, simparams); 
                [P_target, testDV] = calc_covariance_tcmdv(x, t, t_s, stm_t, TCMr_time_test, vel_disp_flag, deltaV, P_i, range, simparams); 
    
                % Comparing
                if testDV < minDV
                    if i == length(TCMr_time_best) && sqrt(trace(P_target(1:3,1:3))) > simparams.P_max_r
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

    

    if isempty(find(gradient_vector~=0))
        improving = 0;
    else
        % Taking a step in the gradient direction
        TCMr_idx_best = TCMr_idx_best + gradient_vector;
        TCMr_time_best = t(TCMr_idx_best)';


        % Storing a new min DV after the step
%         [~, minDV] = calc_covariance_tcmdv_v2(x, t, t_s, stm_t, stm_t_i, TCMr_time_best, vel_disp_flag, deltaV, P_i, range, simparams); 
        [~, minDV] = calc_covariance_tcmdv(x, t, t_s, stm_t, TCMr_time_best, vel_disp_flag, deltaV, P_i, range, simparams); 
        TCMr_idx_test_save(end,:) = TCMr_idx_best;
        minDV_save(end) = minDV;

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
%             outline = [counter, t(TCMr_idx_test_save(minIdx,:))'*ndTime2hrs, minDV*ndVel2kms*1000];
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
        TCMr_time_best = t(TCMr_idx_best)';


        counter = counter + 1;

    end


end

% if output
%     fclose(fileID);
%     ppp=1;
% end


end