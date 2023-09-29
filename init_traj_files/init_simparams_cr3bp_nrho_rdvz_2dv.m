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
simparams.sig_tcm_error = .01 / 1e3 / ndDist2km * ndTime2sec; % Velocity 1 sigma = 1 cm/s
% simparams.sig_tcm_error = .1 / 1e3 / ndDist2km * ndTime2sec; % Velocity 1 sigma = 10 cm/s

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
simparams.Qt = 1e-8 * eye(3);
% simparams.Qt = zeros(3,3);

%% Load saved trajectory parameters

%% Trajectory parameter structure
simparams.m = 7; % number of elements per trajectory segment (6 element state vector, 1 for time duration of segment)
simparams.n = 8; % number of trajectory segments
% simparams.n = 10; % number of trajectory segments

simparams.x0 = zeros(simparams.m, simparams.n); % empty storage for initial trajectory guess

%% Trajectory options

% simparams.maneuverSegments = [2, simparams.n+1]; % the segments with an impulsive maneuver at their beginning (or the nodes with an impulsive maneuver)
simparams.maneuverSegments = [2, simparams.n]; % the segments with an impulsive maneuver at their beginning (or the nodes with an impulsive maneuver)
simparams.P_constrained_nodes = simparams.maneuverSegments(2:end); % Nodes where the position dispersion is constrained to simparams.P_max_r
simparams.max_num_TCMs = 5; % maximum number of TCMs per TCM optimization portion (between nominal maneuvers)

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

simparams.correct_nominal_dvs = 0; % flag to incorporate a dispersion correction with the nominal delta Vs or not

simparams.perform_correction = 1; % flag to incorporate TCM in the trajectory or not

simparams.constrain_dv1_inclination_change = 0; % flag to constrain all inclination change to happen at dv1

simparams.rdvz_flag = 0; % flag to identify a flexible final state for the rendezvous problem setup 

simparams.start_P_growth_node = maneuverSegments(1); % At which node to allow the covariance to grow via linear dynamics/STM. Another way to think about it: where simparams.P_initial is applied initially

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

[~, x0_nrho] = singleDifferentialCorrection(x0_nrho, designMap, constraintMap, T_nrho/2, 0, simparams);






%% Propagate the NRHO to get chaser and target initial states
% time_past_perilune_chaser0 = .25 / ndTime2days;
time_past_perilune_chaser0 = 2.25 / ndTime2days;
% time_past_perilune_chaser0 = .5 / ndTime2days;
% time_past_perilune_chaser0 = 1.5 / ndTime2days;
% time_past_chaser_target0 = 1 / ndTime2hrs; 
time_past_chaser_target0 = 1 * 5/60 / ndTime2hrs; 

[x0_chaser] = stateProp(x0_nrho, time_past_perilune_chaser0, simparams);
[x0_target] = stateProp(x0_chaser, time_past_chaser_target0, simparams);

simparams.x0_target = x0_target;
simparams.T_target = T_nrho;

% Adding these, but they seem superfluous
simparams.x_init = x0_chaser;
simparams.T0 = T_nrho;



% Initial transfer duration
% total_transfer_duration = .7475 / ndTime2days;

% total_transfer_duration = 1.84 / ndTime2days;
total_transfer_duration = 1.8 / ndTime2days;

% total_transfer_duration = 5.6 / ndTime2days;
% total_transfer_duration = 5.3 / ndTime2days;
% total_transfer_duration = 4 / ndTime2days;
simparams.tf = total_transfer_duration;


% % Use a different transfer duration to generate the initial guess
% total_transfer_duration = 4 / ndTime2days;


% Propagate target by total_transfer_duration
simparams.x_target = stateProp(x0_target, total_transfer_duration, simparams);

% Propagate chaser by total_transfer_duration
[xf_chaser, x_nrho_t, t_nrho] = stateProp(x0_chaser, total_transfer_duration, simparams);




% Extra coast in the target orbit duration
% extra_target_coast = .5 / ndTime2days;
% extra_target_coast = .1 / ndTime2days;

extra_target_coast = .25 / ndTime2days;
% extra_target_coast = .21 / ndTime2days;

% extra_target_coast = 2.1 / ndTime2days;



%% Use differential correction to calculate change in velocity for the
% chaser to intersect X_velocity_final

% Free variables - velocity of chaser at t0
designEl = logical([0 0 0 1 1 1]);

% Constraint vector - r_target at t_final minus r_chaser at t_final
constraintEl = logical([1 1 1 0 0 0]);

% Differential correction targeting



X_curr = x0_chaser;
T = total_transfer_duration;

tDesignFlag = 0;
correction = 1;
while correction

    [x_final, stm, xstm_t] = stateStmProp(X_curr, total_transfer_duration, simparams);






    % Construct state vector with STM appended
%     X = [X_curr; reshape(stm_initial, 36, 1)]; 

    % Propagate state and STM for delta_t to rendezvous
%     [~,X] = ode45(@cr3bp_de_s, [0,t_rndvz], X, options);

    % Extract initial state vector at t=0
%     X_initial = X(1,1:6)';
    % Extract final state vector at t=t_final
%     X_final = X(end,1:6)';
    % Extract STM from t=0 to t=t_final
%     stm = reshape(X(end,7:42),6,6);
    
    
    % Design vector - the state elements at t=0 that are being modified to 
    % try and satisfy constraints at t=T
    X_design = X_curr(designEl);
    
    % If time / Period is a design variable, append T to the end of the
    % design vector
    if tDesignFlag == 1
        X_design = [X_design; T];
        Tindex = length(X_design);
    end    
    
    % Constraint vector - the elements being minimized 
    fX_constraint = x_final(constraintEl) - simparams.x_target(constraintEl);

    
    % DF matrix - portions of the STM
    % Partial of constraints wrt design variables  
    % Said another way - the sensitivity of the final constraint variables
    % to variations in the initial design variables    
    DF = stm(constraintEl,designEl);
    
    % If time is a design variable, append the time derivatives of the
    % constraint variables to the DF matrix. 
    % Incorporates the sensitivity of the constraint variables to time
    if tDesignFlag == 1
        % Calculate derivatives wrt time of state at t=T
%         dXdt = cr3bp_sFrame_nd_stm_de(x_final);
        dXdt = stateDot(x_final, mu, simparams.dynSys);
%         dXdt = dXdt(1:6);
        % Append to DF matrix
        DF = [DF, dXdt(constraintEl)];
    end
    
    % Performing corrections to minimize fX_constraint
    % X_next is the adjusted design elements (corr. to initial state, time
    % if included; more options exist)
    X_next = X_design - DF' * inv(DF*DF') * fX_constraint;

    % Re-assign X_curr and T after correction
%     X_curr = zeros(6,1);
%     X_curr(designEl) = X_next(1:length(designEl));
    X_curr(designEl) = X_next(1:sum(designEl));
    
    if tDesignFlag == 1
        T = X_next(Tindex);
    end

    % Check if fX_constraint is sufficiently small
    disp(strcat('Norm of constraint vector after correction iteration:',num2str(norm(fX_constraint))))
    if norm(fX_constraint) < 1e-10        
        X_return = X_curr;
        correction = 0;
    end


end

% X_chaser_mod_0 = X_return;





[xf_chaser_dv1, x_chaser_postdv1, t_chaser_postdv1] = stateProp(X_return, total_transfer_duration, simparams);






num_rdvz_segments = simparams.maneuverSegments(2) - simparams.maneuverSegments(1);



simparams.x0(1:6,1) = x0_chaser;
simparams.x0(7,1) = 1e-10;

x_chaser_coast_postdv1 = subdivide_segment(x_chaser_postdv1, t_chaser_postdv1, num_rdvz_segments);
simparams.x0(:,simparams.maneuverSegments(1):simparams.maneuverSegments(2)-1) = x_chaser_coast_postdv1;


% 
if simparams.maneuverSegments(end) == simparams.n
    % Then need to add coast to the target state
    x_target_coast_start = simparams.x_target;
    simparams.x_target = stateProp(x_target_coast_start, extra_target_coast, simparams);
    simparams.tf = simparams.tf + extra_target_coast;

    simparams.x0(:,end) = [x_target_coast_start; extra_target_coast];



elseif simparams.maneuverSegments(end) == simparams.n + 1
    % Then don't need a coast to the target state

end




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

