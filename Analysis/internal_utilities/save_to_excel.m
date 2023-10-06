% Use the writetable function
data_out = table;
data_out.savename = savename;
if ~isfield(data_out, 'Scenario')
    scenario = 'NRHO rendezvous';
end
data_out.Scenario = scenario;
data_out.dir = outputPath;
data_out.Robust = simparams.perform_correction;

deltaVmps = deltaV*ndVel2kms*1000;
tcmtotalmps = 3*min_tcm_dv*ndVel2kms*1000;
totalcostmps = totalDV*ndVel2kms*1000;

data_out.TotalDV_mps = totalcostmps;
data_out.Nominal_DV_mps = deltaVmps;
data_out.TCM_cost_mps = tcmtotalmps;
data_out.Number_of_segments = simparams.n;
data_out.Initial_position_dispersion_per_axis_km = simparams.sig_pos * ndDist2km;
data_out.Initial_velocity_dispersion_per_axis_mps = simparams.sig_vel * ndDist2km / ndTime2sec * 1000;
data_out.Max_target_position_dispersion_rss = simparams.P_max_r * ndDist2km;
data_out.TCM_execution_error_mps = simparams.sig_tcm_error * ndDist2km / ndTime2sec * 1000;
data_out.Nominal_DV_execution_error_mps = simparams.sig_dv_error * ndDist2km / ndTime2sec * 1000;
data_out.Q_per_axis_m2_s3 = simparams.Qt(1,1) / ndTime2sec^3 * ndDist2m^2;
data_out.Correcting_nominal_maneuvers = simparams.correct_nominal_dvs;
data_out.Dynamical_system = simparams.dynSys;
writetable(data_out,'./sims/robust_results_summary.xlsx','WriteMode','append','WriteRowNames',false);
% writetable(data_out,'./sims/robust_results_summary.xlsx','WriteMode','append');