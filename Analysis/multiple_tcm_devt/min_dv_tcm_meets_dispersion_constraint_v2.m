function [tcm_feasibleMin_time,tcm_feasibleMin_idx,minTcmDV_meets_constraint] = min_dv_tcm_meets_dispersion_constraint_v2(x, t, t_s, stm_t, vel_disp_flag, deltaV, P_i, range, simparams)
%min_dv_tcm_meets_dispersion_constraint Finds all feasible TCM options that
%meet the target dispersion constraint, then selects the lowest DV option
%of those and returns the time and time index.



% tcm_dv_t = zeros(1,length(t));
% tcm_dv_t = zeros(1, range(2) - range(1) + 1);
% tcm_dv_t(end) = 1e8; % the loop below skips the final time index...a correction once at the final time doesn't make sense.
% rP_tcm_time_t = zeros(1,length(t));
% rP_tcm_time_t = zeros(1, range(2) - range(1) + 1);


% % Loop through the time options, save the target position dispersion and
% % the TCM delta V RSS if execute at that time
% % for i = 1:length(t) - 1
% for i = 1:range(2) - range(1)
%     tcm_idx = range(1) + i -1;
%     tcm_time_i = t(tcm_idx);
% %     tcm_time_i = t(i);
% 
%     [Pn_tcm_time_i, tcm_dv_t(i)] = calc_covariance_tcmdv_v2(x, t, t_s, stm_t, stm_t_i, tcm_time_i, vel_disp_flag, deltaV, P_i, range, simparams);
%     rP_tcm_time_t(i) = sqrt(trace(Pn_tcm_time_i(1:3,1:3)));
% end

% Start from the end instead:

iter = 1;
% tcm_idx = range(2)-range(1);
tcm_idx = length(t)-1;
stepSize = 25;
i = 1;
mode = 1; % Mode 1 = searching in reverse from the end in larger increments for the position dispersion constraint to be violated
while iter
    tcm_time(i) = t(tcm_idx);
%     [Pn_tcm_time_i, tcm_dv_t(i)] = calc_covariance_tcmdv_v2(x, t, t_s, stm_t, stm_t_i, tcm_time(i), vel_disp_flag, deltaV, P_i, range, simparams);
    [Pn_tcm_time_i, tcm_dv_t(i)] = calc_covariance_tcmdv(x, t, t_s, stm_t, tcm_time(i), vel_disp_flag, deltaV, P_i, range, simparams);
    rP_tcm_time_t(i) = sqrt(trace(Pn_tcm_time_i(1:3,1:3)));
    
    

    if mode == 1
        if rP_tcm_time_t(i) > simparams.P_max_r
            if i == 1
                % If this is already violated, it won't cross this
                % threshold going in reverse. Probably in this case because
                % something is wrong.
                assert(0);

            else
                mode = 2; % Mode 2 = searching forward, counting by ones, from the point where the dispersion constraint was violated until it is no longer violated
            end
        else   
        
            tcm_idx = tcm_idx - stepSize;
            if tcm_idx < 1
                tcm_idx = 1;
                mode = 2; % a TCM at the first index satisfies the constraint
            end


        
            i = i+1;
        end
    elseif mode == 2
        if rP_tcm_time_t(i) <= simparams.P_max_r
            iter = 0;
        else
            tcm_idx = tcm_idx + 1;
            i = i+1;
        end
    end
end


tcm_feasibleMin_time = t(tcm_idx);

tcm_feasibleMin_idx = tcm_idx;
minTcmDV_meets_constraint = tcm_dv_t(end);

% 
% % Of all the final dispersion options, find the ones that meet the target
% % position dispersion constraint (1 km for here)
% 
% rP_meets_constraint_logical = rP_tcm_time_t <= simparams.P_max_r;
% 
% rP_meets_constraint_logical(end) = false;
% 
% % Find the one that is the minimum delta V
% % first, find the deltaV options that meet the constraint
% tcm_total_t_meets_constraint = tcm_dv_t(rP_meets_constraint_logical);
% 
% 
% % DEBUG IF STATEMENT - if it isn't meeting the constraint, pick a really
% % expensive TCM toward the end of the trajectory...although, faulty logic
% % currently because the final TCM time is fixed in subsequent steps. ACTION
% % 
% if isempty(tcm_total_t_meets_constraint)
%     tcm_feasibleMin_idx = length(t) - 3;
%     tcm_feasibleMin_time = t(tcm_feasibleMin_idx);
%     
%     minTcmDV_meets_constraint = tcm_dv_t(tcm_feasibleMin_idx);
% 
% 
%     assert(0); % until I come back and keep working on this
%     
% else
% 
% 
%     
%     
%     % the find the minimum, and the index of the minimum
%     [minTcmDV_meets_constraint, min_meets_idx] = min(tcm_total_t_meets_constraint);
%     % find the corresponding time index of all tcm options that meet the constraint
%     min_meets_idx = min_meets_idx + range(1) - 1;
%     t_meets = t(rP_meets_constraint_logical);
%     % get the time from the feasible set index and feasible set time array
%     tcm_feasibleMin_time = t_meets(min_meets_idx);
%     tcm_feasibleMin_idx = find(t == tcm_feasibleMin_time);
% 
% 
% end % end DEBUG if statement

end

% 
% % COME BACK HERE AND SHOW DR GELLER FIGURE 3
% 
% %%%%%% DEBUGGING:
% % Covariance growth without any maneuvers or modifications
% P_t = tmult(stm_t, tmult(simparams.P_initial, stm_t, [0 1]));
% 
% for i = 1:size(P_t,3)
%     Pmag_t(i) = sqrt(trace(P_t(1:3,1:3,i)));
% end
% 
% % Covariance growth with a nominal maneuver
% P_t_dvError(:,:,1) = simparams.P_initial;
% 
% for i = 1:length(t(range(1):range(2)))
%     if i > 1
%         P_t_dvError(:,:,i) = calc_covariance_tcmdv(x, t(1:i), t_s(1:i), stm_t(:,:,1:i), [], vel_disp_flag, P_i, simparams);
%     end
%     Pmag_t_dvError(i) = sqrt(trace(P_t_dvError(1:3,1:3,i)));
% end
% 
% 
% % Pmag_t_dvError(end)*ndDist2km
% % Pmag_t(end)*ndDist2km
% 
% figure;
% plot(t(range(1):range(2)),Pmag_t(range(1):range(2)))
% hold on
% plot(t(range(1):range(2)),Pmag_t_dvError(range(1):range(2)))
% xlim([t(range(1)) t(range(2))]);
% 
% title('Position dispersion through traj without TCM')
% 
% 
% figure
% plot(t(range(1):range(2)-1),tcm_dv_t(range(1):range(2)-1))
% % ylim([0 tcm_dv_t(end-100)*1.5])
% xlim([t(range(1)) t(range(2))]);
% title('TCM magnitude as a function of execution time')
% 
% 
% figure
% plot(t(range(1):range(2)-1),rP_tcm_time_t(range(1):range(2)-1))
% xlim([t(range(1)) t(range(2)-1)]);
% hold on
% yline(simparams.P_max_r)
% plot(t(range(1):range(2)-1),rP_meets_constraint_logical(range(1):range(2)-1) * simparams.P_max_r,'LineWidth',2)
% title('Target position dispersion as a function of TCM execution time with constraint shown')
% 
% figure
% 
% plot(t(rP_meets_constraint_logical),tcm_dv_t(rP_meets_constraint_logical),'LineWidth',3)
% xlim([t(1) t(end)]);
% hold on;
% plot(tcm_feasibleMin_time, minTcmDV_meets_constraint,'.','MarkerSize',25)
% title('TCM magnitude options that meet the target position dispersion constraint')
% 
% 
% 
% % Adding a second TCM
% 
% for i = 1:tcm_feasibleMin_idx - 1
%     tcm2_time_i = sort([tcm_feasibleMin_time, t(i)]);
%     if length(tcm2_time_i) == length(unique(tcm2_time_i))
%         [Pn_2tcm_time_i, tcm2_dv_t(i)] = calc_covariance_tcmdv(x, t, t_s, stm_t, tcm2_time_i, vel_disp_flag, P_i, simparams);
%         rP_2tcm_time_t(i) = sqrt(trace(Pn_2tcm_time_i(1:3,1:3)));
%     end
% end
% 
% 
% figure
% plot(t(1:tcm_feasibleMin_idx - 1), tcm2_dv_t)
% title('Second TCM magnitude options as a function of execution time along traj')
% 
% 
% 
% 
% 
% 
% figure
% plot(t(1:tcm_feasibleMin_idx - 1), rP_2tcm_time_t)
% hold on
% yline(simparams.P_max_r)
% 
% 
% 
% 
% 
% 
% 
% 
% % Plot this portion of the trajectory, something weird is going on
% [~,x_t] = stateProp(x(1:6), t(end) - t(1), simparams);
% figure
% axis equal
% plot3(x_t(:,1), x_t(:,2), x_t(:,3))


