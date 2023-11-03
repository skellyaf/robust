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

% Dynamical system flag for which differential equations to use
simparams.dynSys = 'cr3bp'; % options are 2bp and cr3bp.

% Numerical integration options
% simparams.options = odeset('AbsTol',1e-13,'RelTol',1e-13);
simparams.options = odeset('AbsTol',2.3e-14,'RelTol',2.3e-14);


%% Random variables / initial uncertainty
% Units km, km/hr, km/hr^2

% simparams.P_max_r = 100 / ndDist2km; % km converted to ND dist
simparams.P_max_r = 1 / ndDist2km; % km converted to ND dist

% Initial uncertainty
% zero
% simparams.sig_pos = 1e-12;
% simparams.sig_vel = 1e-12;

% Small
% simparams.sig_pos = 10 / 1e3 / ndDist2km; % Position +/- 10 m in all 3 direction
% simparams.sig_vel = 10 / 1e6 / ndDist2km * ndTime2sec; % Velocity +/- 1 cm/s in all 3 directions

% Medium
% simparams.sig_pos = 1 / ndDist2km; % Position +/- 1 km in all 3 direction converted to ND dist
% simparams.sig_vel = 1 / 1e3 / ndDist2km * ndTime2sec; % Velocity +/- 1 m/s in all 3 directions converted to ND dist / ND time

% Large
simparams.sig_pos = 10 / ndDist2km; % Position +/- 10 km in all 3 direction converted to ND dist
simparams.sig_vel = 10 / 1e3 / ndDist2km * ndTime2sec; % Velocity +/- 10 m/s in all 3 directions converted to ND dist / ND time
% simparams.sig_vel = 10 / 1e5 / ndDist2km * ndTime2sec; % Velocity +/- 10 cm/s in all 3 directions converted to ND dist / ND time


simparams.P_initial = diag([simparams.sig_pos^2 simparams.sig_pos^2 simparams.sig_pos^2 simparams.sig_vel^2 simparams.sig_vel^2 simparams.sig_vel^2]);

% TCM execution error
% simparams.sig_tcm_error = 1e-12; % Velocity 1 sigma = nearly 0 cm/s
% simparams.sig_tcm_error = .01 / 1e3 / ndDist2km * ndTime2sec; % Velocity 1 sigma = 1 cm/s
simparams.sig_tcm_error = .1 / 1e3 / ndDist2km * ndTime2sec; % Velocity 1 sigma = 10 cm/s

simparams.R = diag([simparams.sig_tcm_error, simparams.sig_tcm_error, simparams.sig_tcm_error]).^2;

% Nominal maneuver execution error
% simparams.sig_dv_error = 1e-12; % Velocity 1 sigma = nearly 0 cm/s
% simparams.sig_dv_error = .1 / 1e3 / ndDist2km * ndTime2sec; % Velocity 1 sigma = 10 cm/s
simparams.sig_dv_error = 10 / 1e3 / ndDist2km * ndTime2sec; % Velocity 1 sigma = 10 m/s


simparams.R_dv = diag([simparams.sig_dv_error, simparams.sig_dv_error, simparams.sig_dv_error]).^2;


simparams.add_tcm_improvement_threshold = sqrt(trace(simparams.R)) * 3;

% simparams.R = diag([0 0 0]);

% simparams.Qt = sqrt(4.8e-7^2 / 3) * eye(3) * (ndTime2sec^3/ndDist2km^2) * .000001; % the value used for dev/testing
% simparams.Qt = sqrt(4.8e-7^2 / 3) * eye(3) * (ndTime2sec^3/ndDist2km^2) * .00001;
% simparams.Qt = 4.8e-7 * eye(3) * ndTime2sec^3 / ndDist2m^2 * 1;
% simparams.Qt = 1e-6 * eye(3);
simparams.Qt = 1e-8 * eye(3);

%% Load saved trajectory parameters

%% Trajectory parameter structure
simparams.m = 7; % number of elements per trajectory segment (6 element state vector, 1 for time duration of segment)
simparams.n = 30; % number of trajectory segments

simparams.x0 = zeros(simparams.m, simparams.n); % empty storage for initial trajectory guess

%% Trajectory options

simparams.maneuverSegments = [2, 13, simparams.n]; % the segments with an impulsive maneuver at their beginning
simparams.P_constrained_nodes = simparams.maneuverSegments(2:end); % Nodes where the position dispersion is constrained to simparams.P_max_r
simparams.max_num_TCMs = 10; % maximum number of TCMs per TCM optimization portion (between nominal maneuvers)

simparams.nom_dvctied = 0; % 1 A flag to force the TCM to occur concurrently with the corresponding nominal impulsive maneuver identified by the following variable
simparams.maneuver_w_corr = 0; % 1 The index of simpar.maneuverSegments where a correction occurs (currently the first nominal maneuver); set to 0 to not tie to a nominal maneuver

%%% THE FOLLOWING NOT CURRENTLY IMPLEMENTED!!!
simparams.fixed_xfer_duration = 0; % Flag for a fixed total transfer duration (total transfer duration stored in .tf)
% simparams.tf is set inside generateInitialXfer.m;

simparams.idv_tcmR_method = 0; % A flag to use either the traditional position covariance trace method vs the unit vector approximation when the tcm and delta V are concurrent
simparams.idv_tcmV_method = 0;

simparams.seg1_coast_fraction = 0.5; % percent of orbital period to coast in the first segment for initial segment guess
%  TODO: add if statement to determine the coast period if there's a
%  maneuver at the beginning of the first segment. same for last.
simparams.segn_coast_fraction = 0.3; % percent of orbital period to coast in the last segment for initial segment guess

% To target zero position dispersion at the final maneuver (beginning of
% final segment if the final segment is a coast, for example), set the
% following flag to anything but zero:
simparams.target_final_maneuver = 1;

simparams.perform_correction = 1; % flag to incorporate TCM in the trajectory or not
simparams.target_Pr_constraint_on = 1; % Flag to constrain the target position dispersion (relevant when the TCMs are tied to nodes instead of optimized each iteration)


simparams.correct_nominal_dvs = 1; % flag to incorporate a dispersion correction with the nominal delta Vs or not

simparams.constrain_dv1_inclination_change = 0; % flag to constrain all inclination change to happen at dv1

simparams.start_P_growth_node = 2; % At which node to allow the covariance to grow via linear dynamics/STM. Another way to think about it: where simparams.P_initial is applied initially

simparams.tcm_rss_factor = 3;



%% Orbit parameters
%% Initial orbit - currently circular inclined
altitude_init = 450; % km
coe_init.a = altitude_init + earth.rad; % semimajor axis, km
coe_init.ecc = 0.0001; % eccentricity
coe_init.inc = .1 * pi/180; % inclination, deg
coe_init.raan = 180 * pi/180; % right ascension, deg (not used for equatorial orbits)
coe_init.argp = 0; % arg of perigee, deg (not used for circular & equatorial orbits)

% From shooting / poincare map test
nu = 1.8;
coe_init.nu = nu; % true anomaly, deg (not used for circular orbits)

% if initial orbit is circular/equatorial, the following are used.
coe_init.truelon = nan; % deg, angle between x-axis and satellite (only for circular and equatorial)
coe_init.arglat = 225; % deg, angle between ascending node and satellite position (only for circular inclined orbits)
coe_init.lonper = nan; % deg, angle between x-axis and eccentricity vector (only for eccentric equatorial orbits)
simparams.coe_init = coe_init;


% From shooting / poincare map test
dv = 3:.1:3.5;
dv_sol = 2.85;

% Create state from best shooting result
[r_2b_leo,v_2b_leo] = orbel2rv(coe_init.a, coe_init.ecc, coe_init.inc, coe_init.raan, coe_init.argp, coe_init.nu, earth.mu);
r_2b_leo_nd = r_2b_leo / ndDist2km;
r_3b_leo = r_2b_leo_nd - [simparams.mu; 0; 0]; % position of earth is (-mu, 0, 0)
v_2b_leo_nd = v_2b_leo / ndVel2kms;
X_depart_leo_pre_mnvr = [r_3b_leo; v_2b_leo_nd];

i_leo_v = X_depart_leo_pre_mnvr(4:6) / norm(X_depart_leo_pre_mnvr(4:6));
v_leo_nd = X_depart_leo_pre_mnvr(4:6) + dv_sol * i_leo_v;
x_3b_leo_depart_nd = [X_depart_leo_pre_mnvr(1:3); v_leo_nd];



%% Adding in more than 1 full orbit coast in LEO with a deltaV1=0
% the LEO departure Delta V is now DV2, lunar orbit arrival is DV3

simparams.T0 = 2*pi*sqrt(simparams.coe_init.a^3/earth.mu) / ndTime2sec;
% Add coast (in reverse) to get to the DV1 state 
x_leoFullOrbit_start = stateProp(X_depart_leo_pre_mnvr, -simparams.T0, simparams);

% Add an additional coast (in reverse) to get the initial state

simparams.x_init = stateProp(x_leoFullOrbit_start, - simparams.T0 * simparams.seg1_coast_fraction, simparams);

% Add the initial state to trajectory parameter vector
simparams.x0(1:6,1) = simparams.x_init;
simparams.x0(7,1) = simparams.T0 * simparams.seg1_coast_fraction;

% Separate the full orbit LEO coast into multiple segments per simparams.maneuverSegments
[~,x_leoFullOrbit_t, t_leoFullOrbit] = stateProp(x_leoFullOrbit_start, simparams.T0, simparams);

n_leo_orbit = simparams.maneuverSegments(2) - simparams.maneuverSegments(1);
x_leoFullOrbit = subdivide_segment(x_leoFullOrbit_t, t_leoFullOrbit, n_leo_orbit);

simparams.x0(:,simparams.maneuverSegments(1):simparams.maneuverSegments(2)-1) = x_leoFullOrbit;


%% Target orbit - low lunar orbit
altitude_targ = 1000; % Lunar altitude, km
coe_targ.a = altitude_targ + moon.rad; % semimajor axis, km
coe_targ.ecc = 0.001; % eccentricity
% Inclination of 27 degrees (27, 50, 76, & 86 enable extended LLO stays)
% coe_targ.inc = (180 - 27) * pi/180; % inclination, deg
coe_targ.inc = pi;
coe_targ.raan = 180 * pi/180; % Right ascension - 180 degrees (not used for equatorial orbits)
coe_targ.argp = 0; % arg of perigee, deg (not used for circular & equatorial orbits)
coe_targ.nu = 180 * pi/180; % true anomaly, deg (not used for circular orbits)
% if initial orbit is circular/equatorial, the following are used.
coe_targ.truelon = 90; % deg, angle between x-axis and satellite (only for circular and equatorial)
coe_targ.arglat = 180; % deg, angle between ascending node and satellite position (only for circular inclined orbits)
coe_targ.lonper = nan; % deg, angle between x-axis and eccentricity vector (only for eccentric equatorial orbits)
simparams.coe_targ = coe_targ;

% 2 body COEs for LLO
[r_2b_llo,v_2b_llo] = orbel2rv(coe_targ.a,coe_targ.ecc,coe_targ.inc,coe_targ.raan,coe_targ.argp,coe_targ.nu,moon.mu);
r_2b_llo_nd = r_2b_llo / ndDist2km;
v_2b_llo_nd = v_2b_llo / ndVel2kms;

% Position of moon in nondimensional CR3BP is (1-mu)
% Velocity - all Y
r_3b_llo_nd = r_2b_llo_nd + [1-simparams.mu; 0; 0];
% CR3BP state at arrival to low lunar orbit
x_3b_llo_arrive_nd = [r_3b_llo_nd; v_2b_llo_nd];

% Add a coast to get to the actual LLO target (non dimensional)
simparams.T_target = 2*pi*sqrt(simparams.coe_targ.a^3/moon.mu) / ndTime2sec;

% Add final segment to trajectory parameter vector
simparams.x0(1:6,end) = x_3b_llo_arrive_nd;
simparams.x0(7,end) = simparams.T_target * simparams.segn_coast_fraction;

% Propagate from arrival to the target
simparams.x_target = stateProp(x_3b_llo_arrive_nd, simparams.T_target * simparams.segn_coast_fraction, simparams);

%% Transfer states / initial guess

simparams.xfer1_duration = 3.25 / ndTime2days;
[~,x_xfer_t, t] = stateProp(x_3b_leo_depart_nd, simparams.xfer1_duration, simparams);
% plot3(x_xfer_t(:,1),x_xfer_t(:,2),x_xfer_t(:,3))

n = simparams.maneuverSegments(3) - simparams.maneuverSegments(2);
X = subdivide_segment(x_xfer_t, t, n);


simparams.x0(:,simparams.maneuverSegments(2):simparams.maneuverSegments(3)-1) = X;

simparams.x0 = simparams.x0(:);


%% Flyby perilune radius
simparams.constrain_flyby_radius = false; % bool, true or false

%% Load prior traj, assign nodes to TCMs, re-assign segment durations
%%

% Load
% load('eed_leo_planar_robust0.mat'); % x_opt
simparams_old = simparams;
% load('eed_leo_planar_robust_tcmNodes1.mat','x_opt','simparams'); % the robust trajectory loaded starting from the file above. re-baselining to reaply # of TCMs / nodes (load includes simparams)
% simparams_old.maneuverSegments = simparams.maneuverSegments;
% simparams_old.P_constrained_nodes = simparams.P_constrained_nodes;
% simparams_old.n = simparams.n;
% simparams = simparams_old;

load('eed_leo_planar_robust0.mat','x_opt'); % the robust trajectory loaded starting from the file above. re-baselining to reaply # of TCMs / nodes (load includes simparams)


simparams.x0 = x_opt;

%% Perform the initial trajectory history propagation

% To calculate the optimal number of TCMs and place a segment at each TCM
% intersection and modify how the optimization problem is constructed.
traj_setup  = createStateStmSttQdQHistory(simparams.x0, simparams);

if isfield(simparams,'rdvz_flag')
    if simparams.rdvz_flag == 1
        x0 = reshape(simparams.x0,simparams.m,simparams.n);
        totalTime = sum(x0(7,:));
        simparams.x_target = stateProp(simparams.x0_target, totalTime, simparams);
    end
end

% Calculate total impulsive delta V for initial guess trajectory
[deltaV0, deltaVs_nom0] = calcDeltaV(simparams.x0, traj_setup.x_i_f, traj_setup.stm_i, simparams);

% Calculate the optimal TCMs
[tcm_time_setup, tcm_idx_setup, min_tcm_dv_setup, ~, ~, tcm_dv_each_setup] = opt_multiple_tcm_wQ(simparams.x0, traj_setup, deltaVs_nom0, simparams); % inputs: x, t, t_s, stm_t, stm_t_i, simparams

simparams.x0 = reshape(simparams.x0,simparams.m,simparams.n);

% figure
% plotMultiSegTraj(simparams.x0, traj0.x_t, traj0.t_s, simparams, tcm_idx0);

%% Re-assign segment durations and the nodes that are tied to TCMs

x_new = zeros(simparams.m, length(tcm_time_setup)+2 + 2);
x_new(:,1) = simparams.x0(:,1);

x_new(:,end) = simparams.x0(:,end);

x_new(1:6,2) = simparams.x0(1:6,2);

t_dv2 = sum(simparams.x0(7,1:simparams.P_constrained_nodes(1)-1));
tcm_idxs_p1 = find(tcm_time_setup < t_dv2);
tcm_idxs_p2 = find(tcm_time_setup > t_dv2);


x_new(1:6,3 + length(tcm_idxs_p1)) = simparams.x0(1:6,simparams.P_constrained_nodes(1));


% TCM portion 1
for i = 1:length(tcm_idxs_p1)
    % TCM states
    x_tcm_curr = traj_setup.x_t( tcm_idx_setup(tcm_idxs_p1(i)),1:6)';
    x_new(1:6,i+2) = x_tcm_curr;

    % Segment duration
    x_new(7,i+1) = traj_setup.t(tcm_idx_setup(tcm_idxs_p1(i))) - sum(x_new(7,1:i));

end
t_flyby = sum(simparams.x0(7,1:simparams.P_constrained_nodes(1)-1));
x_new(7, length(tcm_idxs_p1)+2) = t_flyby - tcm_time_setup(length(tcm_idxs_p1));

% TCM portion 2
for i = 1:length(tcm_idxs_p2)
    x_new(1:6,i+length(tcm_idxs_p1) + 3) = traj_setup.x_t( tcm_idx_setup(tcm_idxs_p2(i)),1:6)';


    % Duration
    x_new(7,i+length(tcm_idxs_p1) + 2) = traj_setup.t(tcm_idx_setup(tcm_idxs_p2(i))) - sum(x_new(7,1:i+length(tcm_idxs_p1) + 1));
end

t_dv3 = sum(simparams.x0(7,1:simparams.P_constrained_nodes(2)-1));
x_new(7,end-1) = t_dv3 - tcm_time_setup(end);

% Reassign important simparams
simparams.n = size(x_new,2);

simparams.maneuverSegments = [2, 3 + length(tcm_idxs_p1), simparams.n];

simparams.P_constrained_nodes = simparams.maneuverSegments(2:end);

simparams.tcm_nodes = [tcm_idxs_p1+2, tcm_idxs_p2+3];

simparams.x0 = x_new(:);




% Insert additional segments
[x_new, simparams] = insert_additional_segments(simparams.x0, simparams, 3, 3);
simparams.x0 = x_new;

[x_new, simparams] = insert_additional_segments(simparams.x0, simparams, 6, 2);
simparams.x0 = x_new;

[x_new, simparams] = insert_additional_segments(simparams.x0, simparams, 9, 2);
simparams.x0 = x_new;

[x_new, simparams] = insert_additional_segments(simparams.x0, simparams, 12, 1);
simparams.x0 = x_new;

[x_new, simparams] = insert_additional_segments(simparams.x0, simparams, 20, 2);
simparams.x0 = x_new;



[x_new, simparams] = insert_additional_segments(simparams.x0, simparams, 23, 3);
simparams.x0 = x_new;





%% FOR AN UNKONWN REASON, THE SERVER IS NOT CREATING THE SAME INITIAL TRAJECTORY. 
% saving it here, then loading it to run on the server rather than creating
% it via script
% save('init_eed_planar_tcms_moreError_simparams.mat','simparams')
%%

%% Reset numerical propagation options
simparams.options = odeset('AbsTol',2.3e-14,'RelTol',2.3e-14);



%% Fmincon optimization options

% Undefined algorithm

simparams.optoptions = optimoptions('fmincon');

% # Iterations
simparams.optoptions.MaxFunctionEvaluations = 1e6;
simparams.optoptions.MaxIterations = 1e6;

% Algorithm
simparams.optoptions.Algorithm = 'interior-point';
% simparams.optoptions.Algorithm = 'sqp';

% Specifying gradient flags
simparams.optoptions.SpecifyConstraintGradient = true;
simparams.optoptions.SpecifyObjectiveGradient = true;


% Optimality and constraint satisfaction tolerances
simparams.optoptions.OptimalityTolerance = 1e-6;
simparams.optoptions.ConstraintTolerance = 2e-10 * simparams.n * simparams.m;
simparams.optoptions.StepTolerance = 1e-15; 
% simparams.optoptions.FiniteDifferenceStepSize = 1e-5;

% To use parallel processing
% simparams.optoptions.UseParallel = true;

% Fmincon interior point feasibility mode
% simparams.optoptions.EnableFeasibilityMode = true;
% simparams.optoptions.SubproblemAlgorithm = 'cg';

simparams.optoptions.Display='iter';



% To have matlab check the analytical gradient, if being used
% simparams.optoptions = optimoptions('fmincon','Algorithm','interior-point','MaxFunctionEvaluations',3e5,'MaxIterations',1e4,'SpecifyObjectiveGradient',true,'CheckGradients',true,'FiniteDifferenceType','central','FiniteDifferenceStepSize',1e-8);

