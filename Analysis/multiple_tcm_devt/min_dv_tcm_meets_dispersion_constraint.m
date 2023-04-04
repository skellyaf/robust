function [tcm_feasibleMin_time,tcm_feasibleMin_idx,minTcmDV_meets_constraint] = min_dv_tcm_meets_dispersion_constraint(t, stm_t, vel_disp_flag, P_i, simparams)
%min_dv_tcm_meets_dispersion_constraint Finds all feasible TCM options that
%meet the target dispersion constraint, then selects the lowest DV option
%of those and returns the time and time index.

% % OPTIONAL INPUT: P_i
% if length(varargin) < 1
%     P_i = simparams.P_initial;
% else
%     P_i = varargin{1};
% end


tcm_dv_t = zeros(1,length(t));
tcm_dv_t(end) = 1e8; % the loop below skips the final time index...a correction once at the final time doesn't make sense.
rP_tcm_time_t = zeros(1,length(t));


% Loop through the time options, save the target position dispersion and
% the TCM delta V RSS if execute at that time
for i = 1:length(t) - 1
    tcm_time_i = t(i);
    [Pn_tcm_time_i, tcm_dv_t(i)] = calc_covariance_tcmdv(t, stm_t, tcm_time_i, vel_disp_flag, P_i, simparams);
    rP_tcm_time_t(i) = sqrt(trace(Pn_tcm_time_i(1:3,1:3)));
end

% Of all the final dispersion options, find the ones that meet the target
% position dispersion constraint (1 km for here)

rP_meets_constraint_logical = rP_tcm_time_t <= simparams.P_max_r;

rP_meets_constraint_logical(end) = false;

% Find the one that is the minimum delta V
% first, find the deltaV options that meet the constraint
tcm_total_t_meets_constraint = tcm_dv_t(rP_meets_constraint_logical);
% the find the minimum, and the index of the minimum
[minTcmDV_meets_constraint, min_meets_idx] = min(tcm_total_t_meets_constraint);
% find the corresponding time index of all tcm options that meet the constraint
t_meets = t(rP_meets_constraint_logical);
% get the time from the feasible set index and feasible set time array
tcm_feasibleMin_time = t_meets(min_meets_idx);
tcm_feasibleMin_idx = find(t == tcm_feasibleMin_time);




end



%%%%%% DEBUGGING:
% 
% P_t = tmult(stm_t, tmult(simparams.P_initial, stm_t, [0 1]));
% 
% for i = 1:size(P_t,3)
%     Pmag_t(i) = sqrt(trace(P_t(1:3,1:3,i)));
% end
% figure;
% plot(t,Pmag_t)
% xlim([t(1) t(end)]);
% 
% 
% figure
% plot(t,tcm_dv_t)
% % ylim([0 tcm_dv_t(end-100)*1.5])
% xlim([t(1) t(end)]);
% 
% 
% figure
% plot(t,rP_tcm_time_t)
% xlim([t(1) t(end)]);
% hold on
% yline(simparams.P_max_r)
% plot(t,rP_meets_constraint_logical * simparams.P_max_r,'LineWidth',2)
% 
% 
% figure
% 
% plot(t(rP_meets_constraint_logical),tcm_dv_t(rP_meets_constraint_logical),'LineWidth',3)
% xlim([t(1) t(end)]);
% hold on;
% plot(tcm_feasibleMin_time, minTcmDV_meets_constraint,'.','MarkerSize',25)
% 
% 


