
% Load some saved data
load('colorblind_colormap.mat')
load('orbital_params.mat');

% Colorblind colormap
simparams.colorblind = colorblind;

%% BR4BP preamble - Sun-B1 rotating frame

moon.mu = 4902.8;

% Mean motion of the earth-moon system
n_em = sqrt((earth.mu + moon.mu)/moon.a^3);

% Mean motion of the sun-B1 system (characteristic time)
n_sb1 = sqrt((earth.mu + moon.mu + sun.mu)/earth.a^3);
simparams.n_sb1 = n_sb1;

% Orbital period in days = 2*pi/n_sb1 / 24/3600
% Characteristic distance
l_star = earth.a;
simparams.l_star = l_star;
% Characteristic mass
m_star = earth.mass + moon.mass + sun.mass;
simparams.m_star = m_star;

% Conversions from dim to ND - Sun-B1 Frame!
ndTime2sec = 1/n_sb1;
ndTime2hrs = 1/n_sb1/3600;
ndTime2days = 1/n_sb1/3600/24;
ndDist2km = l_star;
ndDist2m = l_star * 1000;
ndVel2kms = l_star * n_sb1;
ndVel2ms = l_star * n_sb1 * 1000;

% mu_underbar
simparams.mub = (earth.mass + moon.mass) / m_star;
mub = simparams.mub;

% mu
simparams.mu = moon.mass/(earth.mass + moon.mass);
mu = simparams.mu;

% BR4BP parameters
simparams.a4 = earth.a / moon.a;
a4 = simparams.a4;

% Dynamical system flag for which differential equations to use
simparams.dynSys = 'br4bp_sb1'; 

% Numerical integration options
simparams.options = odeset('AbsTol',1e-12,'RelTol',1e-12);
% simparams.options = odeset('AbsTol',2.3e-14,'RelTol',2.3e-14);

% Equals the same as n_se = sqrt((earth.mu + moon.mu + sun.mu)/earth.a^3); n_se / n 
% Earth orbital period in days = 2*pi/n_se * ndTime2days

% Earth-Moon line rotation rate in the Sun-B1 rotating frame
simparams.theta_em_dot = -1 + n_em*ndTime2sec;

%% State properties
simparams.m = 8; % one more state variable - theta_em
simparams.nsv = simparams.m - 1; % number of state variables
% simparams.n = 10;

%% Random variables / initial uncertainty
% Units km, km/hr, km/hr^2

% simparams.P_max_r = 100 / ndDist2km; % km converted to ND dist
% simparams.P_max_r = 1 / ndDist2km; % km converted to ND dist
% simparams.P_max_r = .5 / ndDist2km; % km converted to ND dist
simparams.P_max_r = 1 / ndDist2km; % km converted to ND dist
% simparams.P_max_r = .1 / ndDist2km; % km converted to ND dist
% simparams.P_max_r = 5 / ndDist2km; % km converted to ND dist

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

simparams.P_initial = zeros(7,7);
simparams.P_initial(1:6,1:6) = diag([simparams.sig_pos^2 simparams.sig_pos^2 simparams.sig_pos^2 simparams.sig_vel^2 simparams.sig_vel^2 simparams.sig_vel^2]);

% TCM execution error
% simparams.sig_tcm_error = 1e-12; % Velocity 1 sigma = nearly 0 cm/s
% simparams.sig_tcm_error = .01 / 1e3 / ndDist2km * ndTime2sec; % Velocity 1 sigma = 1 cm/s
simparams.sig_tcm_error = .1 / 1e3 / ndDist2km * ndTime2sec; % Velocity 1 sigma = 10 cm/s

simparams.R = diag([simparams.sig_tcm_error, simparams.sig_tcm_error, simparams.sig_tcm_error]).^2;

% Nominal maneuver execution error
% simparams.sig_dv_error = 1e-12; % Velocity 1 sigma = nearly 0 cm/s
% simparams.sig_dv_error = .1 / 1e3 / ndDist2km * ndTime2sec; % Velocity 1 sigma = 10 cm/s
% simparams.sig_dv_error = 1 / 1e3 / ndDist2km * ndTime2sec; % Velocity 1 sigma = 1 m/s
simparams.sig_dv_error = 10 / 1e3 / ndDist2km * ndTime2sec; % Velocity 1 sigma = 10 m/s
% simparams.sig_dv_error = 30 / 1e3 / ndDist2km * ndTime2sec; % Velocity 1 sigma = X m/s

simparams.R_dv = diag([simparams.sig_dv_error, simparams.sig_dv_error, simparams.sig_dv_error]).^2;

simparams.add_tcm_improvement_threshold = sqrt(trace(simparams.R)) * 3;

% simparams.R = diag([0 0 0]);



% simparams.Qt = sqrt(4.8e-7^2 / 3) * eye(3) * (ndTime2sec^3/ndDist2km^2) * .000001; % the value used for dev/testing
% simparams.Qt = sqrt(4.8e-7^2 / 3) * eye(3) * (ndTime2sec^3/ndDist2km^2) * .00001;
% simparams.Qt = 4.8e-7 * eye(3) * ndTime2sec^3 / ndDist2m^2 * 1;
% simparams.Qt = 1e-9 * eye(3);
simparams.Qt = 1e-8 * eye(3);
% simparams.Qt = 1e-6 * eye(3);
% simparams.Qt = zeros(3,3);





%% Orbit parameters

simparams.constrain_circ_x0xf = 0;
simparams.fixed_initial_radius = (100 + earth.rad) / ndDist2km;
simparams.fixed_final_radius = (100 + moon.rad) / ndDist2km;

simparams.mu_earth_nd =  earth.mu / ndDist2km^3 * ndTime2sec^2;
simparams.mu_moom_nd = moon.mu / ndDist2km^3 * ndTime2sec^2;


%% Flyby perilune radius
simparams.constrain_flyby_radius = false; % bool, true or false
% simparams.flyby_radius = (moon.rad + 100) / ndDist2km; % distance from the center of the moon
% simparams.flyby_node = simparams.maneuverSegments(2); % the node that is constrained to a certain distance from the moon

%% Setup
% Initial earth-moon angle
theta_em0 = 24 * pi/180;
 
xe = 1 - mub - 1/a4 * mu * cos(theta_em0);
ye = -1/a4 * mu *sin(theta_em0);

xm = 1 - mub + 1/a4 * (1-mu) * cos(theta_em0);
ym = 1/a4 * (1-mu) * sin(theta_em0);
r_earth = [xe; ye; 0];
r_b1 = [1-mub; 0; 0];


%% State propagation

% dt = 77.75 / ndTime2days;
dt = 80 / ndTime2days;

%% Create initial BLT state

nu_earth = -30*pi/180; % initial earth orbit true anomaly
altitude = 100;

% thrust = 0.108246666666667;
thrust = .108247037037037;
[r_sun_sc, v_sc_sb1_nd] = calc_v_sc_earthOrbit_br4bp(earth, nu_earth, altitude, theta_em0, simparams);
v_sc_depart = v_sc_sb1_nd + [-sin(nu_earth) * thrust; cos(nu_earth) * thrust; 0];


% Initial BLT state
x2 = [r_sun_sc; v_sc_depart; theta_em0];
% Propagate
[x_f, stm, xstm_t, ti] = stateStmProp(x2, dt, simparams);

% Create transfer segment parameters
s2 = [x2; dt];

%% Propagate backwards from x1 to get x0 and simparams.x_init

% dt_backwards_from_1 = .1 / ndTime2hrs;
dt_backwards_from_1 = 0;

% Remove the delta V from x1
x1_nodv = [r_sun_sc; v_sc_sb1_nd; theta_em0];
[x1, ~, xstm_t_0] = stateStmProp(x1_nodv, -dt_backwards_from_1, simparams);

s1 = [x1; dt_backwards_from_1];

%% Add initial lunar orbit that comes close to connecting

theta_emf = x_f(7);
[r_sun_sc_lunarOrbit, v_sc_sb1_lunarOrbit] = calc_v_sc_lunarOrbit_br4bp(moon, 270*pi/180, 100, theta_emf, simparams);

% Propagate lunar orbit
xtarg_lunar = [r_sun_sc_lunarOrbit; v_sc_sb1_lunarOrbit; theta_emf];
% dt_lunar = .5 / ndTime2hrs;
dt_lunar = 0 / ndTime2hrs;
[xf_lunar, ~, xstm_t_lunar] = stateStmProp(xtarg_lunar, dt_lunar, simparams);


simparams.x_target = xstm_t_lunar(end,1:simparams.nsv)';
s3 = [xtarg_lunar; dt_lunar];

simparams.x_init = x1;

%% Assemble trajectory
simparams.x0 = [s1; s2; s3];
simparams.n = 3;



%% Trajectory options

% simparams.maneuverSegments = [2, simparams.n+1]; % the segments with an impulsive maneuver at their beginning (or the nodes with an impulsive maneuver)
simparams.maneuverSegments = [2, simparams.n]; % the segments with an impulsive maneuver at their beginning (or the nodes with an impulsive maneuver)
simparams.P_constrained_nodes = simparams.maneuverSegments(2:end); % Nodes where the position dispersion is constrained to simparams.P_max_r
simparams.max_num_TCMs = 10; % maximum number of TCMs per TCM optimization portion (between nominal maneuvers)

simparams.nom_dvctied = 0; % 1 A flag to force the TCM to occur concurrently with the corresponding nominal impulsive maneuver identified by the following variable
simparams.maneuver_w_corr = 0; % 1 The index of simpar.maneuverSegments where a correction occurs (currently the first nominal maneuver); set to 0 to not tie to a nominal maneuver

simparams.fixed_xfer_duration = 0; % Flag for a fixed total transfer duration (total transfer duration stored in .tf)

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
simparams.perform_correction = 0; % flag to incorporate TCM in the trajectory or not
simparams.target_Pr_constraint_on = 0; % Flag to constrain the target position dispersion (relevant when the TCMs are tied to nodes instead of optimized each iteration)


simparams.correct_nominal_dvs = 1; % flag to incorporate a dispersion correction with the nominal delta Vs or not


simparams.constrain_dv1_inclination_change = 0; % flag to constrain all inclination change to happen at dv1

simparams.rdvz_flag = 0; % flag to identify a flexible final state for the rendezvous problem setup 

% simparams.start_P_growth_node = simparams.maneuverSegments(1); % At which node to allow the covariance to grow via linear dynamics/STM. Another way to think about it: where simparams.P_initial is applied initially
simparams.start_P_growth_node = 1; % At which node to allow the covariance to grow via linear dynamics/STM. Another way to think about it: where simparams.P_initial is applied initially

simparams.num_time_idxs_add_per_seg = 0; % in some cases, the indices around the optimal TCMs are somewhat far apart. Adding some to see if it helps.

simparams.tcm_rss_factor = 3;

%% Orbit parameters

simparams.constrain_circ_x0xf = 1;
simparams.fixed_initial_radius = (100 + earth.rad) / ndDist2km;
simparams.fixed_final_radius = (100 + moon.rad) / ndDist2km;

simparams.mu_earth_nd =  earth.mu / ndDist2km^3 * ndTime2sec^2;
simparams.mu_moom_nd = moon.mu / ndDist2km^3 * ndTime2sec^2;

%% Insert additional segments
[x_new, simparams] = insert_additional_segments(simparams.x0, simparams, 2, 12);
[x_new, simparams] = insert_additional_segments(x_new, simparams, simparams.n-1, 3);
[x_new, simparams] = insert_additional_segments(x_new, simparams, simparams.n-1, 2);
simparams.x0 = x_new;



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
simparams.optoptions.ConstraintTolerance = 1e-10 * simparams.n * simparams.m;
simparams.optoptions.StepTolerance = 1e-15; 
% simparams.optoptions.FiniteDifferenceStepSize = 1e-5;

% To use parallel processing
% simparams.optoptions.UseParallel = true;

% Fmincon interior point feasibility mode
% simparams.optoptions.EnableFeasibilityMode = true;
% simparams.optoptions.SubproblemAlgorithm = 'cg';

simparams.optoptions.Display='iter';

