function [TCMr_idx_best, fevals] = random_unit_search_wQ(x, traj, TCMr_idx_best, vel_disp_flag, deltaV, P_i, minDV, simparams)
%random_unit_search Creates a random unit vector to modify the time index
%elements of the TCM index array to search for an improved delta V solution


% % OPTIONAL INPUT: P_i
% if length(varargin) < 1
%     P_i = simparams.P_initial;
% else
%     P_i = varargin{1};
% end




maneuver_seg = simparams.maneuverSegments;
maneuver_idxs = nan*ones(1,length(maneuver_seg));
for i = 1:length(maneuver_seg)
    match = find(traj.t_s==maneuver_seg(i),1)-1;
    if ~isempty(match)
        maneuver_idxs(i) = match;
    end
end







startIdx = 1;

fevals = 0;

TCMr_time_best = traj.t(TCMr_idx_best)';
% [~, minDV] = calc_covariance_tcmdv_v2(x, t, t_s, stm_t, stm_t_i, TCMr_time_best, vel_disp_flag, deltaV, P_i, range, simparams); 

mod_logical = logical(zeros(1,length(TCMr_idx_best)));
mod_logical(1:end) = true; % Does modify the final TCM index
% mod_logical(1:end-1) = true; % Does not modify the final TCM index

lengthMod = sum(mod_logical);


% Adding on a randomized modification method
improving = 1;
notImprovedCount = 0;

while improving

    TCMr_idx_test = TCMr_idx_best;


    % Method 1 - mod all moddable elements with random vector
%     r1 = randi([-1 1],1,lengthMod);
%     TCMr_idx_test(mod_logical) = TCMr_idx_best(mod_logical) + r1;

    % Method 2 - mod only 2 elements at a time, selected randomly
    %%%% edit - to modify by either 1 or 2
%     r1 = randi([0 1], 1, 2);
%     r1(r1==0) = -1;


    elements_to_mod = randperm(lengthMod,2)+startIdx-1;



    select_from = setdiff(-2:2,0);
    r1 = select_from(randi(length(select_from),2,1));

    mod_logical = logical(zeros(1,length(TCMr_idx_best)));
    mod_logical(elements_to_mod) = true;





    % Modify the vector
    TCMr_idx_test(mod_logical) = TCMr_idx_best(mod_logical) + r1;


    if TCMr_idx_test(1) < 1
        TCMr_idx_test(1) = 1;
    end
    if TCMr_idx_test(end) > length(traj.t)
        TCMr_idx_test(end) = length(traj.t);
    end


    for i = 1:length(maneuver_idxs)
        testloc = find(TCMr_idx_test==maneuver_idxs(i));
        if ~isempty(testloc)
            TCMr_idx_test(testloc) = TCMr_idx_best(testloc);
            
        end
    end

    % Check / remove any indices that are less than 1
    lt1_idx = TCMr_idx_test<1;

    TCMr_idx_test(lt1_idx) = 1:sum(lt1_idx);
    TCMr_idx_test = sort(TCMr_idx_test);

    % Check / remove any indices that are more than the traj length
    gtL_idx = TCMr_idx_test >= length(traj.t);
    TCMr_idx_test(gtL_idx) = length(traj.t)-1;



    TCMr_time_test = traj.t(TCMr_idx_test)';

    % Check for duplicates
    if length(TCMr_time_test) ~= length(unique(TCMr_time_test))
        ppp=1;
%         x
%         TCMr_time_test
        notImprovedCount = notImprovedCount + 1;
    else
        
        [Pf, testDV] = calc_covariance_wQ_tcmdv(x, traj, TCMr_time_test, vel_disp_flag, deltaV, P_i, simparams);
        fevals = fevals + 1;
        
        if testDV <= minDV && sqrt(trace(Pf(1:3,1:3))) <= simparams.P_max_r
            TCMr_time_best = TCMr_time_test;
            TCMr_idx_best = TCMr_idx_test;
            minDV = testDV;
        else
            notImprovedCount = notImprovedCount + 1;
        end

    end

%     [~, testDV] = calc_covariance_tcmdv_v2(x, t, t_s, stm_t, stm_t_i, TCMr_time_test, vel_disp_flag, deltaV, P_i, range, simparams); 



%     if abs(minDV - 111.537939828886) < .0000001
%         fevals
%         ppp=1;
%     end
    
    if notImprovedCount > 4*lengthMod^3

        improving = 0;
    end


end




end