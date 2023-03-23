%% Script to initilize simulation parameters (simparams)

% Load some saved data
load('colorblind_colormap.mat')
load('orbital_params.mat');

% Colorblind colormap
simparams.colorblind = colorblind;

% mu
simparams.mu = earth.mu * 3600 ^ 2; % units: km^3/hr^2

% Dynamical system flag for which differential equations to use
simparams.dynSys = '2bp'; % options are 2bp and cr3bp.

% Numerical integration options
simparams.options = odeset('AbsTol',1e-12,'RelTol',1e-12);

% Any useful conversions required
% Not currently being used, but could put a conversion to normalized dist.
% km2normdist = 1/(alt_target + earth.rad);

%% Random variables / initial uncertainty
% Units km, km/hr, km/hr^2
simparams.P_max_r = 500/3; % km
% simparams.P_max_r = 125; % km

% Initial uncertainty
% zero
% simparams.sig_pos = 0;
% simparams.sig_vel = 0;
% Small
% simparams.sig_pos = 10 / 1e3; % Position +/- 10 m in all 3 direction
% simparams.sig_vel = .01 / 1e3; % Velocity +/- 1 cm/s in all 3 directions
% Medium
% simparams.sig_pos = 1000 / 1e3; % Position +/- 1 km in all 3 direction
% simparams.sig_vel = 1 / 1e3; % Velocity +/- 1 m/s in all 3 directions
% Large
simparams.sig_pos = 10000 / 1e3; % Position +/- 10 km in all 3 direction
simparams.sig_vel = 10 / 1e3 * 3600; % Velocity +/- 10 m/s in all 3 directions
% Huge
% simparams.sig_pos = 100000 / 1e3; % Position +/- 100 km in all 3 direction
% simparams.sig_vel = 100 / 1e3; % Velocity +/- 100 m/s in all 3 directions

simparams.P_initial = diag([simparams.sig_pos^2 simparams.sig_pos^2 simparams.sig_pos^2 simparams.sig_vel^2 simparams.sig_vel^2 simparams.sig_vel^2]);

%% Trajectory vector structure
simparams.m = 7; % number of elements per trajectory segment (6 element state vector, 1 for time duration of segment)
simparams.n = 3; % number of trajectory segments
simparams.x0 = zeros(simparams.n * simparams.m,1); % empty storage for initial trajectory guess

%% Trajectory options

simparams.maneuverSegments = [2, simparams.n]; % the segments with an impulsive maneuver at their beginning

simparams.nom_dvctied = 0; % 1 A flag to force the TCM to occur concurrently with the corresponding nominal impulsive maneuver identified by the following variable
simparams.maneuver_w_corr = 0; % 1 The index of simpar.maneuverSegments where a correction occurs (currently the first nominal maneuver); set to 0 to not tie to a nominal maneuver

simparams.idv_tcmR_method = 0; % A flag to use either the traditional position covariance trace method vs the unit vector approximation when the tcm and delta V are concurrent
simparams.idv_tcmV_method = 0;

simparams.fixed_xfer_duration = 0; % Flag for a fixed total transfer duration (total transfer duration stored in .tf)
% simparams.tf is set inside generateInitialXfer.m;


simparams.seg1_coast_fraction = .5; % percent of orbital period to coast in the first segment for initial segment guess
%  TODO: add if statement to determine the coast period if there's a
%  maneuver at the beginning of the first segment. same for last.
simparams.segn_coast_fraction = .3; % percent of orbital period to coast in the last segment for initial segment guess


% To target zero position dispersion at the final maneuver (beginning of
% final segment if the final segment is a coast, for example), set the
% following flag to anything but zero:
simparams.target_final_maneuver = 1;

simparams.perform_correction = 1; % flag to incorporate TCM in the trajectory or not

simparams.constrain_dv1_inclination_change = 0; % flag to constrain all inclination change to happen at dv1


%% Orbit parameters
% Initial orbit - currently circular inclined
altitude_init = 450; % km
coe_init.a = altitude_init + earth.rad; % semimajor axis, km
coe_init.ecc = 0; % eccentricity
coe_init.inc = 63.4; % inclination, deg
coe_init.raan = 0; % right ascension, deg (not used for equatorial orbits)
coe_init.argp = 199; % arg of perigee, deg (not used for circular & equatorial orbits)
coe_init.nu = 199; % true anomaly, deg (not used for circular orbits)
% if initial orbit is circular/equatorial, the following are used.
coe_init.truelon = nan; % deg, angle between x-axis and satellite (only for circular and equatorial)
coe_init.arglat = 50; % deg, angle between ascending node and satellite position (only for circular inclined orbits)
coe_init.lonper = nan; % deg, angle between x-axis and eccentricity vector (only for eccentric equatorial orbits)
simparams.coe_init = coe_init;

% Target orbit - geostationary
altitude_targ = 894; % km
% altitude_targ = 1200; % km
coe_targ.a = altitude_targ + earth.rad; % semimajor axis, km
coe_targ.ecc = 0; % eccentricity
coe_targ.inc = 99; % inclination, deg
coe_targ.raan = 0; % right ascension, deg (not used for equatorial orbits)
coe_targ.argp = 199; % arg of perigee, deg (not used for circular & equatorial orbits)
coe_targ.nu = 199; % true anomaly, deg (not used for circular orbits)
% if initial orbit is circular/equatorial, the following are used.
coe_targ.truelon = 99; % deg, angle between x-axis and satellite (only for circular and equatorial)
coe_targ.arglat = 50; % deg, angle between ascending node and satellite position (only for circular inclined orbits)
coe_targ.lonper = nan; % deg, angle between x-axis and eccentricity vector (only for eccentric equatorial orbits)
simparams.coe_targ = coe_targ;

%% Fmincon optimization options

% Undefined algorithm
% simparams.optoptions = optimoptions('fmincon','MaxFunctionEvaluations',3e5,'MaxIterations',1e4);

% SQP or Interior point algorithms
% simparams.optoptions = optimoptions('fmincon','Algorithm','sqp','MaxFunctionEvaluations',3e5,'MaxIterations',1e4);
simparams.optoptions = optimoptions('fmincon','Algorithm','interior-point','MaxFunctionEvaluations',3e5,'MaxIterations',1e4);

% Objective function gradient
% simparams.optoptions = optimoptions('fmincon','Algorithm','interior-point','MaxFunctionEvaluations',3e5,'MaxIterations',1e4,'SpecifyObjectiveGradient',true);

% Constraint function gradient
% simparams.optoptions = optimoptions('fmincon','Algorithm','interior-point','MaxFunctionEvaluations',3e5,'MaxIterations',1e4,'SpecifyConstraintGradient',true);
% simparams.optoptions = optimoptions('fmincon','Algorithm','sqp','MaxFunctionEvaluations',3e5,'MaxIterations',1e4,'SpecifyConstraintGradient',true);

%%%%%%%%%%%%%%%%%%% Objective and constraint function gradients 
% simparams.optoptions = optimoptions('fmincon','Algorithm','interior-point','MaxFunctionEvaluations',3e5,'MaxIterations',1e4,'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);
% simparams.optoptions = optimoptions('fmincon','Algorithm','sqp','MaxFunctionEvaluations',3e5,'MaxIterations',1e4,'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);
% outputGradients = 1;
% outputCGradients = 1;

% Optimality and constraint satisfaction tolerances
simparams.optoptions.OptimalityTolerance = 1e-9;
simparams.optoptions.ConstraintTolerance = 1e-7;
% simparams.optoptions.StepTolerance = 1e-10; % use with sqp
% simparams.optoptions.FiniteDifferenceStepSize = 1e-5;

% To use parallel processing
simparams.optoptions.UseParallel = true;

% To have matlab check the analytical gradient, if being used
% simparams.optoptions = optimoptions('fmincon','Algorithm','interior-point','MaxFunctionEvaluations',3e5,'MaxIterations',1e4,'SpecifyObjectiveGradient',true,'CheckGradients',true,'FiniteDifferenceType','central','FiniteDifferenceStepSize',1e-8);

