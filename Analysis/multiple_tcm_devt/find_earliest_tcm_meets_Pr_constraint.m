function [final_tcm_time,final_tcm_idx] = find_earliest_tcm_meets_Pr_constraint(x, t, stm_t ,simparams)
%find_earliest_tcm_meets_Pr_constraint Finds the earliest TCM along that
%trajectory that is required to meet the target position dispersion
%RSS constraint. It is for the purpose of initializing a multiple TCM
%trajectory scenario.


notfound = 1;
i = length(t);

while notfound
    % backwards pass from end of traj
    tcm_time_i = t(i);
    P_tcm_time_i = calc_covariance_history(x, t, stm_t, tcm_time_i, simparams);
    rP_tcm_time_i = sqrt(trace(P_tcm_time_i(1:3,1:3)));

    % initialize comparisons
    if i == length(t)
        last_tcm_time = tcm_time_i;
        last_tcm_idx = i;
    end

    % compare to constraint. if it breaks the constraint, save it
    if rP_tcm_time_i > simparams.P_max_r
        final_tcm_time = last_tcm_time;
        final_tcm_idx = last_tcm_idx;
        notfound = 0;
    end

    last_tcm_time = tcm_time_i;
    last_tcm_idx = i;
    i = i-1;

    % If the max target position constraint is never violated, return -1
    if i == 0
        final_tcm_time = last_tcm_time;
        final_tcm_idx = last_tcm_idx;
        notfound = 0;

    end
end



end