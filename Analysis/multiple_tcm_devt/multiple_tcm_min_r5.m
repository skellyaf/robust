% r3 - r2 ended up incorporating up to 5 TCMs and optimizing the location
% of each via a gradient-like modification to the time index, but each was
% modified individually. attempting to create a gradient vector and modify
% numerous indices simultaneously in r3.


% r2 - functionalized a lot of useful parts.
% working on incorporating a third TCM.
% idea: calc next min starting from first TCM along traj. maybe include the
% first or last TCM in the calc, or both, as fixed time options.

%% For now, plot TCM history and corresponding target position dispersion
% RSS throughout the trajectory
x = simparams.x0;

[stm_i, x_i_f, x_t, stm_t, t, t_s] = createStateStmHistory(x, simparams);

% Using saved dynamics, calculate total impulsive delta V
[deltaV, deltaVs_nom] = calcDeltaV(x, x_i_f, simparams);

% Using saved dynamics, calculate the optimal number of TCMs
tic
[tcm_3sigma,tcm_time, dvR3sigma_tr, dvV3sigma_tr] = tcmPair_rv(x, t, stm_t, deltaVs_nom, simparams);
toc

tic
[tcm_time,tcm_idx,min_tcm_dv,tcm_time_cell,tcm_idx_cell,tcm_num_option_DVs] = opt_multiple_tcm(x, t, stm_t, simparams);
toc

ppp=1;