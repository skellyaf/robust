data_out.NRHO_period_days = simparams.T0 * ndTime2days;
data_out.Chaser_initial_state_days_past_perilune = time_past_perilune_chaser0 * ndTime2days;
data_out.Target_initial_state_hours_past_chaser = time_past_chaser_target0 * ndTime2hrs;
data_out.Total_trajectory_duration_days = simparams.tf * ndTime2days;
x_opt = reshape(x_opt,simparams.m,simparams.n);

data_out.Time_between_deltaVs_days = sum(x_opt(7,simparams.maneuverSegments(1):simparams.maneuverSegments(2)-1)) * ndTime2days;
data_out.Initial_guess_final_coast_days = extra_target_coast * ndTime2days;
writetable(data_out,'./sims/robust_results_nrho_rdvz_summary.xlsx','WriteMode','append','WriteRowNames',false);
% writetable(data_out,'./sims/robust_results_nrho_rdvz_summary.xlsx','WriteMode','append');