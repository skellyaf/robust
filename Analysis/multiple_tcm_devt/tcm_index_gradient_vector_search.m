function [TCMr_idx_best, TCMr_time_best, minDV] = tcm_index_gradient_vector_search(x, t, stm_t, TCMr_idx_best ,startIdx, simparams)
%tcm_index_gradient_search Alters the time index of each individual TCM in
%each direction (earlier and later) to search for an improvement in the
%total delta V

TCMr_time_best = t(TCMr_idx_best)';

improving = 1;

[~, minDV] = calc_covariance_tcmdv(x, t, stm_t, TCMr_time_best, simparams); 

while improving
% for j = 1:nIter

    gradient_vector = zeros(1,length(TCMr_idx_best));
    TCMr_idx_test_save = zeros(length(TCMr_time_best)-1 + 1); % One more than the loop below to add the gradient method to the last element
    minDV_save = ones(1,length(TCMr_time_best)-1 + 1) * 1e8;

    for i = startIdx:length(TCMr_time_best) - 1

        TCMr_time_test = TCMr_time_best;
        TCMr_idx_test = TCMr_idx_best;
        
    
        % Go through each TCM time and modify in each direction (forwards and
        % backwards in time) to see if the delta V total improves. Continue in
        % that direction until it stops improving, update the time in the
        % baseline, then go on to the next.
    
        % Modify the i-th element backwards
%         iter = 1;
        improved = 0;
%         while iter
        if i == 1
            compare = 0;
        else
            compare = TCMr_idx_test(i-1);
        end 

        if TCMr_idx_test(i)-1 > compare
            TCMr_idx_test(i) = TCMr_idx_test(i)-1;
            TCMr_time_test(i) = t(TCMr_idx_test(i));
            [~, testDV] = calc_covariance_tcmdv(x, t, stm_t, TCMr_time_test, simparams); 

            if testDV < minDV
                TCMr_idx_test_save(i,:) = TCMr_idx_test;
                minDV_save(i) = testDV;
%                     TCMr_time_best = TCMr_time_test;
%                     TCMr_idx_best = TCMr_idx_test;
%                     minDV = testDV;

                gradient_vector(i) = -1;

                improved = 1;

%             else
%                 iter = 0; % 
            end
%             end
     
              
        end
        
    
    
    
        % if going backwards didn't result in an improvement, try going
        % forwards instead
    
        if improved == 0
%             iter = 1;
    
%             while iter

            TCMr_idx_test = TCMr_idx_best;
            TCMr_time_test = TCMr_time_best;

            if i == length(gradient_vector)
                compare = 1e8;
            else
                compare = TCMr_idx_test(i+1);
            end

            if TCMr_idx_test(i)+1 < compare
                TCMr_idx_test(i) = TCMr_idx_test(i)+1;
                TCMr_time_test(i) = t(TCMr_idx_test(i));
                [~, testDV] = calc_covariance_tcmdv(x, t, stm_t, TCMr_time_test, simparams); 
    
                if testDV < minDV
                    TCMr_idx_test_save(i,:) = TCMr_idx_test;
                    minDV_save(i) = testDV;
%                     TCMr_time_best = TCMr_time_test;
%                     TCMr_idx_best = TCMr_idx_test;

                    gradient_vector(i) = 1;
                    minDV = testDV;
                else
                    
                    iter = 0; % 
                end
    
            end
    
%             end
    
    
    
        end
    
    
    
    end

    

%     if minDV >= testDV || sum(gradient_vector) == 0
    if isempty(find(gradient_vector~=0))
        improving = 0;
    else
        TCMr_idx_best = TCMr_idx_best + gradient_vector;
        TCMr_time_best = t(TCMr_idx_best)';
        [~, minDV] = calc_covariance_tcmdv(x, t, stm_t, TCMr_time_best, simparams); 
        TCMr_idx_test_save(end,:) = TCMr_idx_best;
        minDV_save(end) = minDV;

        % Check against the saved individual mods for the lowest
        [minDV, minIdx] = min(minDV_save);
        TCMr_idx_best = TCMr_idx_test_save(minIdx,:);
        TCMr_time_best = t(TCMr_idx_best)';

    end


end



end