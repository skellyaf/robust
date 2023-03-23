function [TCMr_time_best,TCMr_idx_best,minDV] = tcm_index_gradient_search(x, t, stm_t, TCMr_time_best, TCMr_idx_best ,nIter, simparams)
%tcm_index_gradient_search Alters the time index of each individual TCM in
%each direction (earlier and later) to search for an improvement in the
%total delta V

[~, minDV] = calc_covariance_tcmdv(x, t, stm_t, TCMr_time_best, simparams); 

for j = 1:nIter
    for i = 2:length(TCMr_time_best)-1

        TCMr_time_test = TCMr_time_best;
        TCMr_idx_test = TCMr_idx_best;
        
    
        % Go through each TCM time and modify in each direction (forwards and
        % backwards in time) to see if the delta V total improves. Continue in
        % that direction until it stops improving, update the time in the
        % baseline, then go on to the next.
    
        % Modify the i-th element backwards
        iter = 1;
        improved = 0;
        while iter
            if TCMr_idx_test(i)-1 > TCMr_idx_test(i-1)
                TCMr_idx_test(i) = TCMr_idx_test(i)-1;
                TCMr_time_test(i) = t(TCMr_idx_test(i)-1);
                [~, testDV] = calc_covariance_tcmdv(x, t, stm_t, TCMr_time_test, simparams); 
    
                if testDV < minDV
                    TCMr_time_best = TCMr_time_test;
                    TCMr_idx_best = TCMr_idx_test;
                    minDV = testDV;
                    improved = 1;
    
                else
                    iter = 0; % 
                end
            end
     
              
        end
        
    
    
    
        % if going backwards didn't result in an improvement, try going
        % forwards instead
    
        if improved == 0
            iter = 1;
    
            while iter
    
                if TCMr_idx_test(i)+1 < TCMr_idx_test(i+1)
                    TCMr_idx_test(i) = TCMr_idx_test(i)+1;
                    TCMr_time_test(i) = t(TCMr_idx_test(i)+1);
                    [~, testDV] = calc_covariance_tcmdv(x, t, stm_t, TCMr_time_test, simparams); 
        
                    if testDV < minDV
                        TCMr_time_best = TCMr_time_test;
                        TCMr_idx_best = TCMr_idx_test;
                        minDV = testDV;
                    else
                        
                        iter = 0; % 
                    end
        
                end
    
            end
    
    
    
        end
    
    
    
    end

end



end