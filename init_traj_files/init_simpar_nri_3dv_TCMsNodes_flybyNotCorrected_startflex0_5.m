%% Script to initialize simulation parameters (simparams)

% Load some saved data
load('colorblind_colormap.mat')
load('orbital_params.mat');

% Colorblind colormap
simparams.colorblind = colorblind;

%% CR3BP preamble

moon.mu = 4902.8;

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

% Earth mu
simparams.mu_earth_nd = earth.mu / ndDist2km^3 * ndTime2sec^2;

% Dynamical system flag for which differential equations to use
simparams.dynSys = 'cr3bp'; % options are 2bp and cr3bp.

% Numerical integration options
% simparams.options = odeset('AbsTol',1e-12,'RelTol',1e-12);
simparams.options = odeset('AbsTol',2.3e-14,'RelTol',2.3e-14);


%% Random variables / initial uncertainty
% Units km, km/hr, km/hr^2

% simparams.P_max_r = 100 / ndDist2km; % km converted to ND dist
simparams.P_max_r = 1 / ndDist2km; % km converted to ND dist
% simparams.P_max_r = 5 / ndDist2km; % km converted to ND dist

% Initial uncertainty
% zero
% simparams.sig_pos = 1e-12;
% simparams.sig_vel = 1e-12;

% Very Small
% simparams.sig_pos = 10 / 1e3 / ndDist2km; % Position +/- 10 m in all 3 direction
% simparams.sig_vel = 10 / 1e6 / ndDist2km * ndTime2sec; % Velocity +/- 1 cm/s in all 3 directions

% Small
simparams.sig_pos = 100 / 1e3 / ndDist2km; % Position +/- 100 m in all 3 direction
% simparams.sig_vel = 10 / 1e6 / ndDist2km * ndTime2sec; % Velocity +/- 1 cm/s in all 3 directions

% Medium
% simparams.sig_pos = 1 / ndDist2km; % Position +/- 1 km in all 3 direction converted to ND dist
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
% simparams.sig_dv_error = .1 / 1e3 / ndDist2km * ndTime2sec; % Velocity 1 sigma = 10 cm/s
simparams.sig_dv_error = 1 / 1e3 / ndDist2km * ndTime2sec; % Velocity 1 sigma = 1 m/s
% simparams.sig_dv_error = 10 / 1e3 / ndDist2km * ndTime2sec; % Velocity 1 sigma = 10 m/s
% simparams.sig_dv_error = 30 / 1e3 / ndDist2km * ndTime2sec; % Velocity 1 sigma = X m/s

simparams.R_dv = diag([simparams.sig_dv_error, simparams.sig_dv_error, simparams.sig_dv_error]).^2;

simparams.add_tcm_improvement_threshold = sqrt(trace(simparams.R)) * 3;

% simparams.R = diag([0 0 0]);


% simparams.Qt = 1e-10 * eye(3) * ndTime2sec^3 / ndDist2m^2; % m^2 / sec^3 converted to ND dist ^ 2 / ND time ^ 3

simparams.Qt = 1e-8 * eye(3) * ndTime2sec^3 / ndDist2m^2; % m^2 / sec^3 converted to ND dist ^ 2 / ND time ^ 3
% simparams.Qt = 1e-6 * eye(3) * ndTime2sec^3 / ndDist2m^2; % m^2 / sec^3 converted to ND dist ^ 2 / ND time ^ 3
% simparams.Qt = (.3e-3)^2 * eye(3) * ndTime2sec^3 / ndDist2m^2; % m^2 / sec^3 converted to ND dist ^ 2 / ND time ^ 3
% simparams.Qt = 1e-5 * eye(3) * ndTime2sec^3 / ndDist2m^2; % m^2 / sec^3 converted to ND dist ^ 2 / ND time ^ 3

%% Load saved trajectory parameters

%% Trajectory parameter structure
simparams.m = 7; % number of elements per trajectory segment (6 element state vector, 1 for time duration of segment)
% simparams.n = 25; % number of trajectory segments
simparams.n = 20; % number of trajectory segments
simparams.nsv = 6; % number of state variables

simparams.x0 = zeros(simparams.m, simparams.n); % empty storage for initial trajectory guess

%% Trajectory options

simparams.maneuverSegments = [2, 10, simparams.n]; % the segments with an impulsive maneuver at their beginning
simparams.P_constrained_nodes = simparams.maneuverSegments(2:end); % Nodes where the position dispersion is constrained to simparams.P_max_r
simparams.corrected_nominal_dvs = logical([1 0 1]); % Logical flag for each nominal maneuver identifying if it should get the combined correction savings
% simparams.correct_nominal_dvs = 0; % flag to incorporate a dispersion correction with the nominal delta Vs or not


simparams.max_num_TCMs = 10; % maximum number of TCMs per TCM optimization portion (between nominal maneuvers)


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

simparams.start_P_growth_node = 2; % At which node to allow the covariance to grow via linear dynamics/STM. Another way to think about it: where simparams.P_initial is applied initially

simparams.tcm_rss_factor = 3;

simparams.target_Pr_constraint_on = 1; % Flag to constrain the target position dispersion (relevant when the TCMs are tied to nodes instead of optimized each iteration)

simparams.skip_dv_1 = false; % Flag to not include the first nominal DV in the cost function / DV calculation / DV gradients.

simparams.circ_init_constraint = true; % flag to constrain the initial orbit to be circular and of a specific energy and radius rather than a full initial state constraint

%% Orbit parameters
%% Initial orbit - currently circular inclined
altitude_init = 450; % km
coe_init.a = altitude_init + earth.rad; % semimajor axis, km
coe_init.ecc = 0.0001; % eccentricity
coe_init.inc = 28 * pi/180; % inclination, deg
coe_init.raan = 50 * pi/180; % right ascension, deg (not used for equatorial orbits)
coe_init.argp = 0; % arg of perigee, deg (not used for circular & equatorial orbits)

% From shooting / poincare map test
coe_init.nu = 3.242; % true anomaly, deg (not used for circular orbits)

% if initial orbit is circular/equatorial, the following are used.
coe_init.truelon = nan; % deg, angle between x-axis and satellite (only for circular and equatorial)
coe_init.arglat = 225; % deg, angle between ascending node and satellite position (only for circular inclined orbits)
coe_init.lonper = nan; % deg, angle between x-axis and eccentricity vector (only for eccentric equatorial orbits)
simparams.coe_init = coe_init;


% From shooting / poincare map test
dv = 2.98;

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
t_tli = 4.0 / ndTime2days;
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

%% Perform single differential correction to ensure it repeats

% Design vector map and elements
% Initial x & z position, y velocity (and appending time T later)
designMap = [1, 0, 1, 0, 1, 0];

% Constraint vector map and elements
% Y position, x velocity, and z velocity all zero at subsequent x-z plane crossing
constraintMap = [0, 1, 0, 1, 0, 1];

[~, x0_nrho] = singleDifferentialCorrection(x0_nrho, designMap, constraintMap, T_nrho/2, 0, simparams);



% Propagate the NRHO
[~,x_nrho_t, t_nrho] = stateProp(x0_nrho, T_nrho*1.015, simparams);

% Indices from external tests
T_coast_nrho_target = t_nrho(280) - t_nrho(230);

[minz,minzidx] = min(x_nrho_t(:,3));

simparams.x_target = x_nrho_t(minzidx,:)';
simparams.T_target = T_nrho;

simparams.x0(1:6,simparams.maneuverSegments(3)) = x_nrho_t(230,:)';
simparams.x0(7,simparams.maneuverSegments(3)) = T_coast_nrho_target;


% 
% 
% 
% figure
% plot3(x_nrho_t(230:280,1),x_nrho_t(230:280,2),x_nrho_t(230:280,3),'LineWidth',2)




% Single parameter vector

% load('nri_det_opt_flex0.mat');
load('nri_det_opt_20seg.mat');
% load('nri_det_opt_update.mat');

simparams.x0 = x_opt;
simparams.maneuverSegments = [2, 10, 20];
simparams.P_constrained_nodes = simparams.maneuverSegments(2:end);
simparams.n = 20;
% simparams.n=size(x_opt,2);
% 
% simparams.maneuverSegments = [2, 14, simparams.n]; % the segments with an impulsive maneuver at their beginning
% simparams.P_constrained_nodes = simparams.maneuverSegments(2:end); % Nodes where the position dispersion is constrained to simparams.P_max_r


%% Perform the initial trajectory history propagation

% To calculate the optimal number of TCMs and place a segment at each TCM
% intersection and modify how the optimization problem is constructed.
traj0  = createStateStmSttQdQHistory(simparams.x0, simparams);

if isfield(simparams,'rdvz_flag')
    if simparams.rdvz_flag == 1
        x0 = reshape(simparams.x0,simparams.m,simparams.n);
        totalTime = sum(x0(7,:));
        simparams.x_target = stateProp(simparams.x0_target, totalTime, simparams);
    end
end

% Calculate total impulsive delta V for initial guess trajectory
[deltaV0, deltaVs_nom0] = calcDeltaV(simparams.x0, traj0.x_i_f, traj0.stm_i, simparams);

% Calculate the optimal TCMs
% [tcm_time0, tcm_idx0, min_tcm_dv0, ~, ~, tcm_dv_each0] = opt_multiple_tcm_wQ(simparams.x0, traj0, deltaVs_nom0, simparams); % inputs: x, t, t_s, stm_t, stm_t_i, simparams
[tcm_time0, tcm_idx0, min_tcm_dv0, ~, ~, tcm_dv_each0] = opt_multiple_tcm_wQ_multiPart(simparams.x0, traj0, deltaVs_nom0, simparams); % inputs: x, t, t_s, stm_t, stm_t_i, simparams

simparams.x0 = reshape(simparams.x0,simparams.m,simparams.n);


%% Re-assign segment distances and the nodes that are tied to TCMs

x_new = zeros(simparams.m, length(tcm_time0)+2 + 2);
x_new(:,1) = simparams.x0(:,1);

x_new(:,end) = simparams.x0(:,end);



% What is the index for the second nominal maneuver / first TCM set target
t_tcm_tgt_1 = sum(simparams.x0(7,1:simparams.P_constrained_nodes(1)-1));
idx_tcm_tgt_1 = find(traj0.t == t_tcm_tgt_1);

% tcms_before_p1 = tcm_idx0 ( tcm_idx0 < idx_tcm_tgt_1)

% tcm_idxs_p1 = 1:4;
% tcm_idxs_p2 = 5:7;
tcm_idxs_p1 = find( tcm_idx0 < idx_tcm_tgt_1);
tcm_idxs_p2 = find( tcm_idx0 > idx_tcm_tgt_1);

x_new(1:6,2) = simparams.x0(1:6,2);
x_new(1:6,length(tcm_idxs_p1) + 3) = simparams.x0(1:6,simparams.P_constrained_nodes(1));

% TCM portion 1
for i = 1:length(tcm_idxs_p1)
    % TCM states
    x_tcm_curr = traj0.x_t( tcm_idx0(tcm_idxs_p1(i)),1:6)';
    x_new(1:6,i+2) = x_tcm_curr;

    % Segment duration
    x_new(7,i+1) = traj0.t(tcm_idx0(tcm_idxs_p1(i))) - sum(x_new(7,1:i));

end
t_flyby = sum(simparams.x0(7,1:simparams.P_constrained_nodes(1)-1));
x_new(7, length(tcm_idxs_p1) + 2) = t_flyby - tcm_time0(length(tcm_idxs_p1));

% TCM portion 2
for i = 1:length(tcm_idxs_p2)
    x_new(1:6,i+length(tcm_idxs_p1) + 3) = traj0.x_t( tcm_idx0(tcm_idxs_p2(i)),1:6)';


    % Duration
    x_new(7,i+length(tcm_idxs_p1) + 2) = traj0.t(tcm_idx0(tcm_idxs_p2(i))) - sum(x_new(7,1:i+length(tcm_idxs_p1) + 1));
end

t_nri = sum(simparams.x0(7,1:simparams.P_constrained_nodes(2)-1));
x_new(7,length(tcm_idx0) + 3) = t_nri - tcm_time0(end);

% Reassign important simparams
simparams.n = size(x_new,2);

simparams.maneuverSegments = [2, 3 + length(tcm_idxs_p1), simparams.n];

simparams.P_constrained_nodes = simparams.maneuverSegments(2:end);

simparams.tcm_nodes = [tcm_idxs_p1+2, tcm_idxs_p2+3];


%% Subdivide any lengthy segments
extend_segs = [];
% extend_segs2 = [];

for i = 2:simparams.n-1
%     if x_new(7,i) > .07
%         extend_segs2 = [extend_segs2, i];
%     end

    if sum(traj0.t_s == i) > 200 || x_new(7,i) > .09
        extend_segs = [extend_segs, i];
    end

end
orig_extend_segs = extend_segs;


% extend_segs = find(x_new(1:end-1)>.07)
x_new_save = x_new;

for i = 1:length(extend_segs)
    x_old = x_new;
    seg_i = extend_segs(i);

    [~,x_i_t, t_i] = stateProp(x_old(1:6,seg_i), x_old(7,seg_i), simparams);
%     [~,~,~,~,~,x_i_t, t_i] = statestmsttQQ2Prop(x_old(1:6,seg_i), x_old(7,seg_i), simparams);

%     n_new_segs = ceil(length(t_i) / 300);
    n_new_segs = ceil(sum(traj0.t_s == orig_extend_segs(i)) / 200);
    n_new_segs2 = ceil(t_i(end) / .05);

    if n_new_segs2 > n_new_segs
        n_new_segs = n_new_segs2;
    end

    x_i_new = subdivide_segment(x_i_t, t_i, n_new_segs);

    x_new = zeros(7,size(x_old,2) + n_new_segs - 1);
    x_new(:,1:seg_i - 1) = x_old(:,1:seg_i - 1);
    x_new(:,seg_i:seg_i + n_new_segs - 1) = x_i_new;
    x_new(:,seg_i + n_new_segs:end) = x_old(:,seg_i + 1:end);






%     t_end_seg_i = sum(x_new(7,1:seg_i + 3));

    simparams.tcm_nodes(simparams.tcm_nodes > seg_i) = simparams.tcm_nodes(simparams.tcm_nodes > seg_i) + n_new_segs - 1;

    extend_segs(extend_segs > seg_i) = extend_segs(extend_segs > seg_i) + n_new_segs - 1;

    simparams.maneuverSegments(simparams.maneuverSegments > seg_i) = simparams.maneuverSegments(simparams.maneuverSegments > seg_i) + n_new_segs - 1;

end

simparams.n = size(x_new,2);
simparams.P_constrained_nodes = simparams.maneuverSegments(2:end);

simparams.x0 = x_new(:);

% % % % % % Seg 4 is long, divide it into 4 segments (+3)
% % % % % [~,x_4_t, t_4] = stateProp(x_new(1:6,4), x_new(7,4), simparams);
% % % % % x_4_new = subdivide_segment(x_4_t, t_4, 4);
% % % % % 
% % % % % x_new2 = zeros(7,size(x_new,2)+3);
% % % % % x_new2(:,1:3) = x_new(:,1:3);
% % % % % x_new2(:,4:7) = x_4_new;
% % % % % x_new2(:,8:end) = x_new(:,5:end);
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % % Re-do params
% % % % % simparams.n = size(x_new2,2);
% % % % % 
% % % % % simparams.maneuverSegments = [simparams.maneuverSegments(1), simparams.maneuverSegments(2:end) + 3];
% % % % % 
% % % % % simparams.P_constrained_nodes = simparams.maneuverSegments(2:end);
% % % % % 
% % % % % simparams.tcm_nodes = [simparams.tcm_nodes(1:2), simparams.tcm_nodes(3:end) + 3];
% % % % % 
% % % % % 
% % % % % 
% % % % % % Seg 3 is also long, divide it into 3 segments (+2)
% % % % % 
% % % % % [~,x_3_t, t_3] = stateProp(x_new(1:6,3), x_new(7,3), simparams);
% % % % % x_3_new = subdivide_segment(x_3_t, t_3, 3);
% % % % % 
% % % % % x_new3 = zeros(7,size(x_new2,2)+2);
% % % % % x_new3(:,1:2) = x_new2(:,1:2);
% % % % % x_new3(:,3:5) = x_3_new;
% % % % % x_new3(:,6:end) = x_new2(:,4:end);
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % % Re-do params
% % % % % simparams.n = size(x_new3,2);
% % % % % 
% % % % % simparams.maneuverSegments = [simparams.maneuverSegments(1), simparams.maneuverSegments(2:end) + 2];
% % % % % 
% % % % % simparams.P_constrained_nodes = simparams.maneuverSegments(2:end);
% % % % % 
% % % % % simparams.tcm_nodes = [simparams.tcm_nodes(1), simparams.tcm_nodes(2:end) + 2];
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % % Seg 15 is also long, divide it into 4 segments (+3)
% % % % % 
% % % % % [~,x_15_t, t_15] = stateProp(x_new3(1:6,15), x_new3(7,15), simparams);
% % % % % x_15_new = subdivide_segment(x_15_t, t_15, 4);
% % % % % 
% % % % % x_new4 = zeros(7,size(x_new3,2)+3);
% % % % % x_new4(:,1:14) = x_new3(:,1:14);
% % % % % x_new4(:,15:18) = x_15_new;
% % % % % x_new4(:,19:end) = x_new3(:,16:end);
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % % Re-do params
% % % % % simparams.n = size(x_new4,2);
% % % % % 
% % % % % simparams.maneuverSegments = [simparams.maneuverSegments(1:2), simparams.maneuverSegments(3:end) + 3];
% % % % % 
% % % % % simparams.P_constrained_nodes = simparams.maneuverSegments(2:end);
% % % % % 
% % % % % simparams.tcm_nodes = [simparams.tcm_nodes(1:7), simparams.tcm_nodes(8:end) + 3];
% % % % % 
% % % % % 
% % % % % simparams.x0 = x_new4(:);



% % Seg 19 is also long, divide it into 3 segments (+2)
% 
% [~,x_19_t, t_19] = stateProp(x_new4(1:6,19), x_new4(7,19), simparams);
% x_19_new = subdivide_segment(x_19_t, t_19, 3);
% 
% x_new5 = zeros(7,size(x_new4,2)+2);
% x_new5(:,1:18) = x_new4(:,1:18);
% x_new5(:,19:21) = x_19_new;
% x_new5(:,22:end) = x_new4(:,20:end);
% 
% 
% 
% 
% % Re-do params
% simparams.n = size(x_new5,2);
% 
% simparams.maneuverSegments = [simparams.maneuverSegments(1:2), simparams.maneuverSegments(3:end) + 2];
% 
% simparams.P_constrained_nodes = simparams.maneuverSegments(2:end);
% 
% simparams.tcm_nodes = [simparams.tcm_nodes(1:8), simparams.tcm_nodes(9:end) + 2];
% 
% 
% 

% simparams.x0 = x_new5(:);



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

