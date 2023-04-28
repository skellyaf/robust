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
ndVel2kms = Rm * n;

% mu
simparams.mu = moon.mu/(earth.mu + moon.mu);
mu = simparams.mu;

% Dynamical system flag for which differential equations to use
simparams.dynSys = 'cr3bp'; % options are 2bp and cr3bp.

% Numerical integration options
simparams.options = odeset('AbsTol',1e-12,'RelTol',1e-12);
% simparams.options = odeset('AbsTol',2.3e-14,'RelTol',2.3e-14);


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

% Large
simparams.sig_pos = 10 / ndDist2km; % Position +/- 10 km in all 3 direction converted to ND dist
simparams.sig_vel = 10 / 1e3 / ndDist2km * ndTime2sec; % Velocity +/- 10 m/s in all 3 directions converted to ND dist / ND time
% simparams.sig_vel = 10 / 1e5 / ndDist2km * ndTime2sec; % Velocity +/- 10 cm/s in all 3 directions converted to ND dist / ND time


simparams.P_initial = diag([simparams.sig_pos^2 simparams.sig_pos^2 simparams.sig_pos^2 simparams.sig_vel^2 simparams.sig_vel^2 simparams.sig_vel^2]);

% TCM execution error
% simparams.sig_tcm_error = 1e-12; % Velocity 1 sigma = nearly 0 cm/s
simparams.sig_tcm_error = .01 / 1e3 / ndDist2km * ndTime2sec; % Velocity 1 sigma = 1 cm/s
% simparams.sig_tcm_error = .1 / 1e3 / ndDist2km * ndTime2sec; % Velocity 1 sigma = 10 cm/s

simparams.R = diag([simparams.sig_tcm_error, simparams.sig_tcm_error, simparams.sig_tcm_error]).^2;

% Nominal maneuver execution error
% simparams.sig_dv_error = 1e-12; % Velocity 1 sigma = nearly 0 cm/s
simparams.sig_dv_error = .1 / 1e3 / ndDist2km * ndTime2sec; % Velocity 1 sigma = 10 cm/s

simparams.R_dv = diag([simparams.sig_dv_error, simparams.sig_dv_error, simparams.sig_dv_error]).^2;

simparams.add_tcm_improvement_threshold = sqrt(trace(simparams.R)) * 3;

% simparams.R = diag([0 0 0]);

%% Load saved trajectory parameters

%% Trajectory parameter structure
simparams.m = 7; % number of elements per trajectory segment (6 element state vector, 1 for time duration of segment)
simparams.n = 25; % number of trajectory segments

simparams.x0 = zeros(simparams.m, simparams.n); % empty storage for initial trajectory guess

%% Trajectory options

simparams.maneuverSegments = [2, 14, simparams.n]; % the segments with an impulsive maneuver at their beginning
simparams.P_constrained_nodes = simparams.maneuverSegments(2:end); % Nodes where the position dispersion is constrained to simparams.P_max_r
simparams.max_num_TCMs = 6; % maximum number of TCMs per TCM optimization portion (between nominal maneuvers)

simparams.nom_dvctied = 0; % 1 A flag to force the TCM to occur concurrently with the corresponding nominal impulsive maneuver identified by the following variable
simparams.maneuver_w_corr = 0; % 1 The index of simpar.maneuverSegments where a correction occurs (currently the first nominal maneuver); set to 0 to not tie to a nominal maneuver

%%% THE FOLLOWING NOT CURRENTLY IMPLEMENTED!!!
simparams.fixed_xfer_duration = 0; % Flag for a fixed total transfer duration (total transfer duration stored in .tf)
% simparams.tf is set inside generateInitialXfer.m;

simparams.idv_tcmR_method = 0; % A flag to use either the traditional position covariance trace method vs the unit vector approximation when the tcm and delta V are concurrent
simparams.idv_tcmV_method = 0;

simparams.seg1_coast_fraction = 0.2; % percent of orbital period to coast in the first segment for initial segment guess
%  TODO: add if statement to determine the coast period if there's a
%  maneuver at the beginning of the first segment. same for last.
simparams.segn_coast_fraction = 0.25; % percent of orbital period to coast in the last segment for initial segment guess

% To target zero position dispersion at the final maneuver (beginning of
% final segment if the final segment is a coast, for example), set the
% following flag to anything but zero:
simparams.target_final_maneuver = 1;

simparams.perform_correction = 1; % flag to incorporate TCM in the trajectory or not

simparams.constrain_dv1_inclination_change = 0; % flag to constrain all inclination change to happen at dv1

%% Orbit parameters
%% Initial orbit - currently circular inclined
altitude_init = 450; % km
coe_init.a = altitude_init + earth.rad; % semimajor axis, km
coe_init.ecc = 0.0001; % eccentricity
coe_init.inc = .1 * pi/180; % inclination, deg
coe_init.raan = 180 * pi/180; % right ascension, deg (not used for equatorial orbits)
coe_init.argp = 0; % arg of perigee, deg (not used for circular & equatorial orbits)

% From shooting / poincare map test
coe_init.nu = .8; % true anomaly, deg (not used for circular orbits)

% if initial orbit is circular/equatorial, the following are used.
coe_init.truelon = nan; % deg, angle between x-axis and satellite (only for circular and equatorial)
coe_init.arglat = 225; % deg, angle between ascending node and satellite position (only for circular inclined orbits)
coe_init.lonper = nan; % deg, angle between x-axis and eccentricity vector (only for eccentric equatorial orbits)
simparams.coe_init = coe_init;


% From shooting / poincare map test
dv = 3;

% Create state from best shooting result
[r_2b_leo,v_2b_leo] = orbel2rv(coe_init.a, coe_init.ecc, coe_init.inc, coe_init.raan, coe_init.argp, coe_init.nu, earth.mu);
r_2b_leo_nd = r_2b_leo / ndDist2km;
r_3b_leo = r_2b_leo_nd - [simparams.mu; 0; 0]; % position of earth is (-mu, 0, 0)
v_2b_leo_nd = v_2b_leo / ndVel2kms;
X_depart_leo_pre_mnvr = [r_3b_leo; v_2b_leo_nd];

i_leo_v = X_depart_leo_pre_mnvr(4:6) / norm(X_depart_leo_pre_mnvr(4:6));
v_leo_nd = X_depart_leo_pre_mnvr(4:6) + dv * i_leo_v;
x_3b_leo_depart_nd = [X_depart_leo_pre_mnvr(1:3); v_leo_nd];




%% Create segments from LEO departure to LLO flyby
% Propagate TLI portion
t_tli = 2.9 / ndTime2days;
[~,x_tli_flyby_t,t_tli_flyby] = stateProp(x_3b_leo_depart_nd, t_tli, simparams);

num_tli_flyby_segments = simparams.maneuverSegments(2) - simparams.maneuverSegments(1);

x_tli_lloFlyby = subdivide_segment(x_tli_flyby_t, t_tli_flyby, num_tli_flyby_segments);

simparams.x0(:,simparams.maneuverSegments(1):simparams.maneuverSegments(2)-1) = x_tli_lloFlyby;



%% Initial coast in LEO
simparams.T0 = 2*pi*sqrt(simparams.coe_init.a^3/earth.mu) / ndTime2sec;

% Add coast (in reverse) to get to the initial state 
simparams.x_init = stateProp(X_depart_leo_pre_mnvr, -simparams.T0 * simparams.seg1_coast_fraction, simparams);

% Add the initial state to trajectory parameter vector
simparams.x0(1:6,1) = simparams.x_init;
simparams.x0(7,1) = simparams.T0 * simparams.seg1_coast_fraction;


%% Generating LLO powered flyby to NRHO trajectory portion
altitude_flyby = 1000; % Lunar altitude, km
coe_flyby.a = altitude_flyby + moon.rad; % semimajor axis, km
coe_flyby.ecc = 0.001; % eccentricity
% Inclination of 27 degrees (27, 50, 76, & 86 enable extended LLO stays)
% coe_targ.inc = (180 - 27) * pi/180; % inclination, deg
coe_flyby.inc = -65 * pi/180;
coe_flyby.raan = 325 * pi/180; % Right ascension - (not used for equatorial orbits)
coe_flyby.argp = 0; % arg of perigee, deg (not used for circular & equatorial orbits)
coe_flyby.nu = 310 * pi/180; % true anomaly, deg (not used for circular orbits)
% if initial orbit is circular/equatorial, the following are used.
coe_flyby.truelon = 90; % deg, angle between x-axis and satellite (only for circular and equatorial)
coe_flyby.arglat = 180; % deg, angle between ascending node and satellite position (only for circular inclined orbits)
coe_flyby.lonper = nan; % deg, angle between x-axis and eccentricity vector (only for eccentric equatorial orbits)
simparams.coe_targ = coe_flyby;

% 2 body COEs for LLO flyby
[r_2b_llo,v_2b_llo] = orbel2rv(coe_flyby.a,coe_flyby.ecc,coe_flyby.inc,coe_flyby.raan,coe_flyby.argp,coe_flyby.nu,moon.mu);
r_2b_llo_nd = r_2b_llo / ndDist2km;
v_2b_llo_nd = v_2b_llo / ndVel2kms;

r_3b_llo = r_2b_llo_nd + [1-mu; 0; 0];
    
X_init = [r_3b_llo; v_2b_llo_nd];
T_init = 1.25/ndTime2hrs;

% Propagate the LLO initial state to find LLO departure state
[~,x_llo_t] = stateProp(X_init, T_init, simparams);

% The LLO departure
j = 10;
x_llo_depart = x_llo_t(j,:)';
T_post_flyby = .72 / ndTime2days;

iv_depart = x_llo_depart(4:6)/norm(x_llo_depart(4:6));
i=23; % from shooting test
dv = i * .025 * iv_depart;

x_flyby_depart = x_llo_depart + [0; 0; 0; dv];

[~,x_flyby_nri_t, t_flyby_nri] = stateProp(x_flyby_depart, T_post_flyby, simparams);

% Divide the posty flyby to NRI into segments
num_flyby_to_nri_segments = simparams.maneuverSegments(3) - simparams.maneuverSegments(2);
x_flyby_to_nri = subdivide_segment(x_flyby_nri_t, t_flyby_nri, num_flyby_to_nri_segments);

% Add to initial traj params
simparams.x0(:,simparams.maneuverSegments(2):simparams.maneuverSegments(3)-1) = x_flyby_to_nri;


%% Flyby perilune radius
simparams.constrain_flyby_radius = true; % bool, true or false
simparams.flyby_radius = (moon.rad + 100) / ndDist2km; % distance from the center of the moon
simparams.flyby_node = simparams.maneuverSegments(2); % the node that is constrained to a certain distance from the moon

%% Add coast in target orbit to traj params
% Get initial state for NRHO
load('lagrangeOrbitICs.mat');
x0_nrho = lagrangeOrbitICs.L2_southern.state_nd(:,1);
T_nrho = lagrangeOrbitICs.L2_southern.T_nd(1);

% Propagate the NRHO
[~,x_nrho_t, t_nrho] = stateProp(x0_nrho, T_nrho*1.015, simparams);

% Indices from external tests
T_coast_nrho_target = t_nrho(280) - t_nrho(230);
simparams.x_target = x_nrho_t(280,:)';
simparams.T_target = T_nrho;

simparams.x0(1:6,simparams.maneuverSegments(3)) = x_nrho_t(230,:)';
simparams.x0(7,simparams.maneuverSegments(3)) = T_coast_nrho_target;


% 
% 
% 
% figure
% plot3(x_nrho_t(230:280,1),x_nrho_t(230:280,2),x_nrho_t(230:280,3),'LineWidth',2)




% Single parameter vector
simparams.x0 = simparams.x0(:);

%% Reset numerical propagation options
simparams.options = odeset('AbsTol',2.3e-14,'RelTol',2.3e-14);


%% Fmincon optimization options

% Undefined optimization algorithm
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
simparams.optoptions.StepTolerance = 1e-14; % use with sqp
% simparams.optoptions.FiniteDifferenceStepSize = 1e-5;

% To use parallel processing
% simparams.optoptions.UseParallel = true;

% Fmincon interior point feasibility mode
% simparams.optoptions.EnableFeasibilityMode = true;
% simparams.optoptions.SubproblemAlgorithm = 'cg';


% To have matlab check the analytical gradient, if being used
% simparams.optoptions = optimoptions('fmincon','Algorithm','interior-point','MaxFunctionEvaluations',3e5,'MaxIterations',1e4,'SpecifyObjectiveGradient',true,'CheckGradients',true,'FiniteDifferenceType','central','FiniteDifferenceStepSize',1e-8);

