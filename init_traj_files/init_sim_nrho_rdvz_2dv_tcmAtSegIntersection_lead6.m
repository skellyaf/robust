%% Script to initialize simulation parameters (simparams)

% Load some saved data
load('colorblind_colormap.mat')
load('orbital_params.mat');

% Colorblind colormap
simparams.colorblind = colorblind;

%% CR3BP preamble


Rm = moon.a; 
n = sqrt((earth.mu + moon.mu)/Rm^3);

% Conversions from dim to ND
ndTime2sec = 1/n;
ndTime2hrs = 1/n/3600;
ndTime2days = 1/n/3600/24;
ndDist2km = Rm;
ndDist2m = Rm * 1000;
ndVel2kms = Rm * n;

% mu
simparams.mu = moon.mu/(earth.mu + moon.mu);
mu = simparams.mu;

% Dynamical system flag for which differential equations to use
simparams.dynSys = 'cr3bp'; % options are 2bp and cr3bp.

% Numerical integration options
% simparams.options = odeset('AbsTol',1e-12,'RelTol',1e-12);
simparams.options = odeset('AbsTol',2.3e-14,'RelTol',2.3e-14);


%% Random variables / initial uncertainty
% Units km, km/hr, km/hr^2

% simparams.P_max_r = 100 / ndDist2km; % km converted to ND dist
% simparams.P_max_r = 1 / ndDist2km; % km converted to ND dist
% simparams.P_max_r = .5 / ndDist2km; % km converted to ND dist
% simparams.P_max_r = 1 / ndDist2km; % km converted to ND dist
simparams.P_max_r = .1 / ndDist2km; % km converted to ND dist
% simparams.P_max_r = 5 / ndDist2km; % km converted to ND dist

% Initial uncertainty
% zero
% simparams.sig_pos = 1e-12;
% simparams.sig_vel = 1e-12;

% Small
% simparams.sig_pos = 10 / 1e3 / ndDist2km; % Position +/- 10 m in all 3 direction
% simparams.sig_vel = 10 / 1e6 / ndDist2km * ndTime2sec; % Velocity +/- 1 cm/s in all 3 directions

% Medium
simparams.sig_pos = 1 / ndDist2km; % Position +/- 1 km in all 3 direction converted to ND dist
simparams.sig_vel = 1 / 1e3 / ndDist2km * ndTime2sec; % Velocity +/- 1 m/s in all 3 directions converted to ND dist / ND time

% Large
% simparams.sig_pos = 10 / ndDist2km; % Position +/- 10 km in all 3 direction converted to ND dist
% simparams.sig_vel = 10 / 1e3 / ndDist2km * ndTime2sec; % Velocity +/- 10 m/s in all 3 directions converted to ND dist / ND time
% simparams.sig_vel = 10 / 1e5 / ndDist2km * ndTime2sec; % Velocity +/- 10 cm/s in all 3 directions converted to ND dist / ND time


simparams.P_initial = diag([simparams.sig_pos^2 simparams.sig_pos^2 simparams.sig_pos^2 simparams.sig_vel^2 simparams.sig_vel^2 simparams.sig_vel^2]);

% TCM execution error
% simparams.sig_tcm_error = 1e-12; % Velocity 1 sigma = nearly 0 cm/s
% simparams.sig_tcm_error = .01 / 1e3 / ndDist2km * ndTime2sec; % Velocity 1 sigma = 1 cm/s
simparams.sig_tcm_error = .1 / 1e3 / ndDist2km * ndTime2sec; % Velocity 1 sigma = 10 cm/s

simparams.R = diag([simparams.sig_tcm_error, simparams.sig_tcm_error, simparams.sig_tcm_error]).^2;

% Nominal maneuver execution error
% simparams.sig_dv_error = 1e-12; % Velocity 1 sigma = nearly 0 cm/s
simparams.sig_dv_error = .1 / 1e3 / ndDist2km * ndTime2sec; % Velocity 1 sigma = 10 cm/s
% simparams.sig_dv_error = 1 / 1e3 / ndDist2km * ndTime2sec; % Velocity 1 sigma = 1 m/s
% simparams.sig_dv_error = 10 / 1e3 / ndDist2km * ndTime2sec; % Velocity 1 sigma = 10 m/s
% simparams.sig_dv_error = 30 / 1e3 / ndDist2km * ndTime2sec; % Velocity 1 sigma = X m/s

simparams.R_dv = diag([simparams.sig_dv_error, simparams.sig_dv_error, simparams.sig_dv_error]).^2;

simparams.add_tcm_improvement_threshold = sqrt(trace(simparams.R)) * 3;

% simparams.R = diag([0 0 0]);



% simparams.Qt = sqrt(4.8e-7^2 / 3) * eye(3) * (ndTime2sec^3/ndDist2km^2) * .000001; % the value used for dev/testing
% simparams.Qt = sqrt(4.8e-7^2 / 3) * eye(3) * (ndTime2sec^3/ndDist2km^2) * .00001;
% simparams.Qt = 4.8e-7 * eye(3) * ndTime2sec^3 / ndDist2m^2 * 1;
% simparams.Qt = 1e-8 * eye(3);
simparams.Qt = 1e-6 * eye(3);
% simparams.Qt = zeros(3,3);

%% Load saved trajectory parameters

%% Trajectory parameter structure
simparams.m = 7; % number of elements per trajectory segment (6 element state vector, 1 for time duration of segment)
simparams.n = 10; % number of trajectory segments
% simparams.n = 10; % number of trajectory segments

simparams.x0 = zeros(simparams.m, simparams.n); % empty storage for initial trajectory guess

%% Trajectory options

% simparams.maneuverSegments = [2, simparams.n+1]; % the segments with an impulsive maneuver at their beginning (or the nodes with an impulsive maneuver)
simparams.maneuverSegments = [2, simparams.n]; % the segments with an impulsive maneuver at their beginning (or the nodes with an impulsive maneuver)
simparams.P_constrained_nodes = simparams.maneuverSegments(2:end); % Nodes where the position dispersion is constrained to simparams.P_max_r
simparams.max_num_TCMs = 3; % maximum number of TCMs per TCM optimization portion (between nominal maneuvers)

simparams.nom_dvctied = 0; % 1 A flag to force the TCM to occur concurrently with the corresponding nominal impulsive maneuver identified by the following variable
simparams.maneuver_w_corr = 0; % 1 The index of simpar.maneuverSegments where a correction occurs (currently the first nominal maneuver); set to 0 to not tie to a nominal maneuver

simparams.fixed_xfer_duration = 1; % Flag for a fixed total transfer duration (total transfer duration stored in .tf)

simparams.idv_tcmR_method = 0; % A flag to use either the traditional position covariance trace method vs the unit vector approximation when the tcm and delta V are concurrent
simparams.idv_tcmV_method = 0;

% simparams.seg1_coast_fraction = 0.2; % percent of orbital period to coast in the first segment for initial segment guess
%  TODO: add if statement to determine the coast period if there's a
%  maneuver at the beginning of the first segment. same for last.
% simparams.segn_coast_fraction = 0.25; % percent of orbital period to coast in the last segment for initial segment guess

% To target zero position dispersion at the final maneuver (beginning of
% final segment if the final segment is a coast, for example), set the
% following flag to anything but zero:
simparams.target_final_maneuver = 1;

simparams.correct_nominal_dvs = 1; % flag to incorporate a dispersion correction with the nominal delta Vs or not

simparams.perform_correction = 1; % flag to incorporate TCM in the trajectory or not

simparams.constrain_dv1_inclination_change = 0; % flag to constrain all inclination change to happen at dv1

simparams.rdvz_flag = 0; % flag to identify a flexible final state for the rendezvous problem setup 

% simparams.start_P_growth_node = simparams.maneuverSegments(1); % At which node to allow the covariance to grow via linear dynamics/STM. Another way to think about it: where simparams.P_initial is applied initially
simparams.start_P_growth_node = 2; % At which node to allow the covariance to grow via linear dynamics/STM. Another way to think about it: where simparams.P_initial is applied initially

simparams.num_time_idxs_add_per_seg = 0; % in some cases, the indices around the optimal TCMs are somewhat far apart. Adding some to see if it helps.

simparams.tcm_rss_factor = 3;

simparams.target_Pr_constraint_on = 1; % Flag to constrain the target position dispersion (relevant when the TCMs are tied to nodes instead of optimized each iteration)



%% Orbit parameters



%% Flyby perilune radius
simparams.constrain_flyby_radius = false; % bool, true or false
% simparams.flyby_radius = (moon.rad + 100) / ndDist2km; % distance from the center of the moon
% simparams.flyby_node = simparams.maneuverSegments(2); % the node that is constrained to a certain distance from the moon

%% Add coast in target orbit to traj params
% Get initial state for NRHO
load('lagrangeOrbitICs.mat');
x0_nrho = lagrangeOrbitICs.L2_southern.state_nd(:,1);
T_nrho = lagrangeOrbitICs.L2_southern.T_nd(1);

%% Perform single differential correction to ensure it repeats

% Design vector map and elements
% Initial x & z position, y velocity (and appending time T later)
designMap = [1, 0, 1, 0, 1, 0];

% Constraint vector map and elements
% Y position, x velocity, and z velocity all zero at subsequent x-z plane crossing
constraintMap = [0, 1, 0, 1, 0, 1];

% x0_nrho(1:3) = x0_nrho(1:3) * ondDist2km;
% x0_nrho(1:3) = x0_nrho(1:3) * ondVel2kms;

[~, x0_nrho] = singleDifferentialCorrection(x0_nrho, designMap, constraintMap, T_nrho/2, 0, simparams);






%% Propagate the NRHO to get chaser and target initial states
time_past_perilune_chaser0 = 1.30464 / ndTime2days;
time_past_chaser_target0 = 1 * 6.06288/60 / ndTime2hrs; 

[x0_chaser] = stateProp(x0_nrho, time_past_perilune_chaser0, simparams);
[x0_target] = stateProp(x0_chaser, time_past_chaser_target0, simparams);

% Extra segment 1 coast
chaser_coast_before_dv1 = .25 / ndTime2days;
[x_chaser_before_dv1] = stateProp(x0_chaser, chaser_coast_before_dv1, simparams);

simparams.x0_target = x0_target;
simparams.T_target = T_nrho;

% Extra coast in the target orbit duration
% extra_target_coast = 1.15 / ndTime2days;
% extra_target_coast = .1 / ndTime2days;
extra_target_coast = .5 / ndTime2days;
% extra_target_coast = .21 / ndTime2days;
% extra_target_coast = 2.1 / ndTime2days;

% Adding these, but they seem superfluous
simparams.x_init = x0_chaser;
simparams.T0 = T_nrho;



% Initial transfer duration
% total_transfer_duration = .7475 / ndTime2days;

% total_transfer_duration = 1.84 / ndTime2days;
% total_transfer_duration = 5.4 / ndTime2days;

% total_transfer_duration = 6. / ndTime2days;
total_transfer_duration = 5.08 / ndTime2days;
% total_transfer_duration = 4 / ndTime2days;
simparams.tf = total_transfer_duration;


% % Use a different transfer duration to generate the initial guess
% total_transfer_duration = 4 / ndTime2days;


% Propagate target by total_transfer_duration minus the extra coast
x_rdvz = stateProp(x0_target, total_transfer_duration - extra_target_coast, simparams);
simparams.x_target = stateProp(x_rdvz, extra_target_coast, simparams);

% Propagate chaser by total_transfer_duration minus the coasts on either
% end
time_between_dvs = total_transfer_duration - chaser_coast_before_dv1 - extra_target_coast;

% [xf_chaser, x_nrho_t, t_nrho] = stateProp(x_chaser_before_dv1, time_between_dvs, simparams);






%% Load previous case 1 3 TCM robust trajectory
% reloading and restarting because the optimal switched to 2 TCMs
% load('nrho_rdvz_robust_3tcm_c2.mat'); % contains x_opt
load('nrho_rdvz_robust_2tcm_c1.mat'); % case 1, contains x_opt
simparams.x0=x_opt;


%% Build a nominal 3 segment trajectory

num_rdvz_segments = simparams.maneuverSegments(2) - simparams.maneuverSegments(1);


% 
% simparams.x0(1:6,1) = x0_chaser;
% simparams.x0(7,1) = chaser_coast_before_dv1;
% 
% x_chaser_coast_postdv1 = subdivide_segment(x_chaser_postdv1, t_chaser_postdv1, num_rdvz_segments);
% simparams.x0(:,simparams.maneuverSegments(1):simparams.maneuverSegments(2)-1) = x_chaser_coast_postdv1;
% simparams.x0(1:6,end) = x_rdvz;
% simparams.x0(7,end) = extra_target_coast;
% 
% 
% % Single parameter vector
% simparams.x0 = simparams.x0(:);

%% Reset numerical propagation options
simparams.options = odeset('AbsTol',2.3e-14,'RelTol',2.3e-14);
% simparams.options = odeset('AbsTol',1e-12,'RelTol',1e-12);
% simparams.options.MaxStep = 60 / ndTime2sec;

%% Perform the initial trajectory history propagations
% To calculate the optimal number of TCMs and place a segment at each TCM
% intersection and modify how the optimization problem is constructed.
simparams.n = length(simparams.x0(:))/simparams.m;
simparams.maneuverSegments = [2, simparams.n];
simparams.P_constrained_nodes = simparams.maneuverSegments(2:end);
traj0  = createStateStmSttQdQHistory(simparams.x0, simparams);

if isfield(simparams,'rdvz_flag')
    if simparams.rdvz_flag == 1
        x0 = reshape(simparams.x0,simparams.m,simparams.n);
        totalTime = sum(x0(7,:));
        simparams.x_target = stateProp(simparams.x0_target, totalTime, simparams);
    end
end

% Calculate total impulsive delta V for initial guess trajectory
[deltaV0, deltaVs_nom0, deletelatergradient] = calcDeltaV(simparams.x0, traj0.x_i_f, traj0.stm_i, simparams);

% Calculate the optimal TCMs
[tcm_time0, tcm_idx0, min_tcm_dv0, ~, ~, tcm_dv_each0] = opt_multiple_tcm_wQ(simparams.x0, traj0, deltaVs_nom0, simparams); % inputs: x, t, t_s, stm_t, stm_t_i, simparams

simparams.x0 = reshape(simparams.x0,simparams.m,simparams.n);

%% Rebuild the trajectory segments for simparams.x0
% Calculate the number of segments: 
% 2 (initial and final coasts) + numTcms + 1 (transfer portion gets divided
% n times into n+1 segments)
simparams.n = 2 + length(tcm_time0) + 1;
x0_new = zeros(simparams.m, simparams.n);
% Assign initial segment
x0_new(:,1) = simparams.x0(:,1);
% Assign final segment
x0_new(:,end) = simparams.x0(:,end);
% Assign from DV1 to TCM 1 segment 
x0_new(1:6,2) = simparams.x0(1:6,2);
x0_new(7,2) = tcm_time0(1) - sum(x0_new(7,1));
% Assign states from TCM 1 to TCM 3
x0_new(1:6,3:3+length(tcm_time0)-1) = traj0.x_t(tcm_idx0,:)';

for i = 1:length(tcm_time0)-1
    x0_new(7,i+2) = tcm_time0(i+1) - sum(x0_new(7,1:i+1));
end

x0_new(7,end-1) = simparams.tf - sum(x0_new(7,:));
simparams.x0 = x0_new(:);


simparams.maneuverSegments = [2, simparams.n]; % the segments with an impulsive maneuver at their beginning (or the nodes with an impulsive maneuver)
simparams.P_constrained_nodes = simparams.maneuverSegments(2:end); % Nodes where the position dispersion is constrained to simparams.P_max_r


simparams.tcm_nodes = [3, 4, 5];





%% Fmincon optimization options

% Undefined optimization algorithm
simparams.optoptions = optimoptions('fmincon');

% # Iterations
simparams.optoptions.MaxFunctionEvaluations = 1e6;
simparams.optoptions.MaxIterations = 1e6;

% Algorithm
simparams.optoptions.Algorithm = 'interior-point';
% simparams.optoptions.Algorithm = 'sqp-legacy';
% simparams.optoptions.Algorithm = 'active-set';

% Specifying gradient flags
simparams.optoptions.SpecifyConstraintGradient = true;
simparams.optoptions.SpecifyObjectiveGradient = true;


% Optimality and constraint satisfaction tolerances
% simparams.optoptions.OptimalityTolerance = 1e-15;
% simparams.optoptions.OptimalityTolerance = 1e-8;
simparams.optoptions.OptimalityTolerance = 1e-6;
% simparams.optoptions.ConstraintTolerance = 2e-10 * simparams.n * simparams.m;
simparams.optoptions.ConstraintTolerance = 1e-14 * simparams.n * simparams.m;
simparams.optoptions.StepTolerance = 1e-16; 
% simparams.optoptions.StepTolerance = 1e-10; 
% simparams.optoptions.FiniteDifferenceStepSize = 1e-5;

% To use parallel processing
% simparams.optoptions.UseParallel = true;

% Fmincon interior point feasibility mode
% simparams.optoptions.EnableFeasibilityMode = true;
% simparams.optoptions.SubproblemAlgorithm = 'cg';

simparams.optoptions.Display='iter';

% To have matlab check the analytical gradient, if being used
% simparams.optoptions = optimoptions('fmincon','Algorithm','interior-point','MaxFunctionEvaluations',3e5,'MaxIterations',1e4,'SpecifyObjectiveGradient',true,'CheckGradients',true,'FiniteDifferenceType','central','FiniteDifferenceStepSize',1e-8);

