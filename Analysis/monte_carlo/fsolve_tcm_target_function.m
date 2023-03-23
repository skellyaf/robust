function [del_r] = fsolve_tcm_target_function(tcm1_dv, simparams, x_opt, deltaVs_nom, maneuver_times, tcm_time, tcm_target_time, x0)

event_times = [maneuver_times,tcm_time, tcm_target_time];
event_dvs = [deltaVs_nom, tcm1_dv, zeros(3,1)];

[~, ~, ~, x_curre] = mcProp(x0,event_times,event_dvs,simparams);


del_xe = x_opt(1:6,simparams.maneuverSegments(end)) - x_curre;

del_r = del_xe(1:3);


end