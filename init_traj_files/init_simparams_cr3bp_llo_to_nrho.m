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
simparams.options = odeset('AbsTol',2.3e-14,'RelTol',2.3e-14);
% simparams.options = odeset('AbsTol',1e-12,'RelTol',1e-12);


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

% simparams.R = diag([0 0 0]);

%% Load saved trajectory parameters

%% Trajectory parameter structure
simparams.m = 7; % number of elements per trajectory segment (6 element state vector, 1 for time duration of segment)
simparams.n = 12; % number of trajectory segments
simparams.x0 = zeros(simparams.m, simparams.n); % empty storage for initial trajectory guess

%% Trajectory options

simparams.maneuverSegments = [2, simparams.n]; % the segments with an impulsive maneuver at their beginning

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






%% Generating low lunar orbit state
a = moon.rad + 110; % Lunar altitude of 110km
e = .001; % Slightly eccentric
% Inclination of 86 degrees (27, 50, 76, & 86 enable extended LLO stays)
inc = 86 * pi/180;
% Right ascension - 180 degrees
W = 180 * pi/180;
% Arg of perigee - zero degrees
w = 0;
% True anomaly - 180 degrees;
nu = 180 * pi/180;

% 2 body COEs for LLO
[r_2b_llo,v_2b_llo] = orbel2rv(a,e,inc,W,w,nu,moon.mu);
r_2b_llo_nd = r_2b_llo / ndDist2km;
v_2b_llo_nd = v_2b_llo / ndVel2kms;

% Position of moon in nondimensional CR3BP is (1-mu)
% Velocity - all Y

r_3b_llo = r_2b_llo_nd + [1-mu; 0; 0];

X_init = [r_3b_llo; v_2b_llo_nd];
T_init = 1.5/ndTime2hrs;

% Test propagate the LLO initial conditions for 2 hours
[~,x_llo_t, t_llo] = stateProp(X_init, T_init, simparams);
simparams.T0 = 2/ndTime2hrs;

% figure
% % plot3(x_llo_t(:,1),x_llo_t(:,2),x_llo_t(:,3),'LineWidth',2)
% xlabel('X (nd)')
% ylabel('Y (nd)')
% zlabel('Z (nd)')
% title('Position in Earth-centered Earth-Moon rotating frame')
% axis equal
% grid on
% hold on
% plot3(X_init(1), X_init(2), X_init(3), '.','MarkerSize',25)


%% Generating NRHO state
addpath('C:\Users\skell\OneDrive - USU\Documents\PhD\research\CR3BP\data');
load('lagrangeOrbitICs.mat');


% I believe this is the 9:2 NRHO, but it has been a while...should double check 
x0_nrho = lagrangeOrbitICs.L2_southern.state_nd(:,1);
T_nrho = lagrangeOrbitICs.L2_southern.T_nd(1);
simparams.T_target = T_nrho;

% Propagate and plot the NRHO
[~,x_nrho_t, t_nrho] = stateProp(x0_nrho, T_nrho*1.015, simparams);

% plot3(x_nrho_t(:,1),x_nrho_t(:,2),x_nrho_t(:,3),'LineWidth',2)







reverse_prop_idx = 400;
x4 = [x0_nrho; t_nrho(reverse_prop_idx)];
% x4 = [x0_nrho; t_nrho(175)];
simparams.x_target = x_nrho_t(reverse_prop_idx,:)';


% plot3(x_nrho_t(400,1),x_nrho_t(400,2),x_nrho_t(400,3),'.','MarkerSize',25)




%% Generating transfer options

x_llo_depart = x_llo_t(end,:)';
T_test = 3 / ndTime2days;

iv_depart = x_llo_depart(4:6)/norm(x_llo_depart(4:6));

% for i = 59:64
%     
%     dv = i * .01 * iv_depart;
% 
%     x_test = x_llo_depart + [0; 0; 0; dv];
%    
%     [~,x_test_t] = stateProp(x_test, T_test, mu,'cr3bp');
% 
%     plot3(x_test_t(:,1),x_test_t(:,2),x_test_t(:,3))
% 
% 
% end

i = 62;
dv = i * .01 * iv_depart;
x_test = x_llo_depart + [0; 0; 0; dv];
[~,x_test_t] = stateProp(x_test, T_test, simparams);
% plot3(x_test_t(:,1),x_test_t(:,2),x_test_t(:,3),'LineWidth',2)


modx = linspace(-.025,.025,5);
mody = linspace(-25,25,10);
modz = linspace(-4,4,10);


for i = 3
    for j = 9
        for k = 1
            dvm = dv;
%             dvm(1) = dv(1);
            dvm(1) = dv(1) + dv(1)*modx(i);
            dvm(2) = dv(2) + dv(2)*mody(j);
            dvm(3) = dv(3) + dv(3)*modz(k);

            x_test = x_llo_t(end-i*2+1,:)' + [0; 0; 0; dvm];
            [~,x_test_t] = stateProp(x_test, T_test, simparams);
%             plot3(x_test_t(:,1),x_test_t(:,2),x_test_t(:,3),'LineWidth',2)
            

            T_init = t_llo(end-i*2+1);

        end
        p=2;
    end
end


x_depart = x_test;

x2 = [x_depart; T_test];

%% Propagate backwards from end state and modify to try and get a connection

xend_nrho = x_nrho_t(end,:)';
T_back = -2 / ndTime2days;


iv_end = -xend_nrho(4:6)/norm(xend_nrho(4:6));

for i = 1
    
    dv = i * .01 * iv_end;

    x_test = xend_nrho + [0; 0; 0; dv];
   
    [~,x_test_t] = stateProp(x_test, T_back, simparams);

%     plot3(x_test_t(:,1),x_test_t(:,2),x_test_t(:,3),'LineWidth',2)


end
% x_depart = x_test;



modx = linspace(-1,1,10);
mody = linspace(-1,1,10);
modz = linspace(-1,1,10);


x_test_back = x_test;

for i = 10
    for j = 5
        for k = 8


            x_test = x_test_back;
            x_test(4) = x_test(4) * (1 + 4*modx(i));
            x_test(5) = x_test(5) * (1 + .1*mody(j));
            x_test(6) = x_test(6) * (1 + .5*modx(k));
            [~,x_test_t] = stateProp(x_test, T_back*.7, simparams);
%             plot3(x_test_t(:,1),x_test_t(:,2),x_test_t(:,3),'LineWidth',2)      

            

        end
        p=2;
    end
end

x3 = [x_test_t(end,:)'; -T_back*.7];


%% try 1 more time to get a better depart matchup


% modx = linspace(-1,1,11);
% mody = linspace(-1,1,11);
% modz = linspace(-1,1,11);

% x_test_fwd = x_depart;

% for i = 6:10
%     for j = 1:11
%         for k = 1:11
% 
% 
%             x_test = x_depart;
%             x_test(4) = x_test(4) * (1 + .1*modx(i));
%             x_test(5) = x_test(5) * (1 + 4*mody(j));
%             x_test(6) = x_test(6) * (1 + .1*modx(k));
%             [~,x_test_t] = stateProp(x_test, T_test*1.5, mu,'cr3bp');
%             plot3(x_test_t(:,1),x_test_t(:,2),x_test_t(:,3),'LineWidth',1)      
% 
%             
% 
%         end
%         p=2;
%     end
%
% end
%


%% Put together state vector
x1 = [X_init; T_init];
simparams.x_init = x1(1:6);
x = [x1; x2; x3; x4];




%% Separate into more segments

n = simparams.maneuverSegments(2) - simparams.maneuverSegments(1);
% Break forward portion up into n/2 segments and reverse portion into n/2

[~,x_2_fwd_t, t_2_fwd] = stateProp(x2(1:6), x2(7), simparams);

p1segs = round(n/2);
p2segs = n - p1segs;

idx_sep = floor(length(t_2_fwd) / p1segs);

iter = 1;
i = idx_sep;
x_2_segs = [];
x_2_segs(1:6,1) = x2(1:6);
x_2_segs(7,1) = t_2_fwd(i);



while iter

    if i > length(t_2_fwd) - idx_sep
        if i - idx_sep > 0
            x_2_segs(7,end) = t_2_fwd(end) - t_2_fwd(i - idx_sep);
        end
        iter = 0;
    else
        x_2_segs(1:6,end+1) = x_2_fwd_t(i,:)';
        x_2_segs(7,end) = t_2_fwd(i + idx_sep) - t_2_fwd(i);
    end
    i = i + idx_sep;   
    



end

% Break reverse portion up into n/2 segments
[~,x_3_t, t_3] = stateProp(x3(1:6), x3(7), simparams);


idx_sep = floor(length(t_3) / p2segs);

iter = 1;
i = idx_sep;
x_3_segs = [];
x_3_segs(1:6,1) = x3(1:6);
x_3_segs(7,1) = t_3(i);



while iter

    if i > length(t_3) - idx_sep
        if i - idx_sep > 0
            x_3_segs(7,end) = t_3(end) - t_3(i - idx_sep);
        end
        iter = 0;
    else
        x_3_segs(1:6,end+1) = x_3_t(i,:)';
        x_3_segs(7,end) = t_3(i + idx_sep) - t_3(i);
    end
    i = i + idx_sep;   
    



end

%% Build X now!
% x0 = [x1, x_2_segs, x_3_segs, x4];
x0 = [x1, x_2_segs, x_3_segs, [simparams.x_target; 0]];
x0(7,11) = .0001;
simparams.x0 = x0(:);






%% Fmincon optimization options

% Undefined algorithm
% simparams.optoptions = optimoptions('fmincon','MaxFunctionEvaluations',3e5,'MaxIterations',1e4);

% SQP or Interior point algorithms
% simparams.optoptions = optimoptions('fmincon','Algorithm','sqp','MaxFunctionEvaluations',3e5,'MaxIterations',1e4);
% simparams.optoptions = optimoptions('fmincon','Algorithm','interior-point','MaxFunctionEvaluations',3e5,'MaxIterations',1e4);

% Objective function gradient
% simparams.optoptions = optimoptions('fmincon','Algorithm','interior-point','MaxFunctionEvaluations',3e5,'MaxIterations',1e4,'SpecifyObjectiveGradient',true);

% Constraint function gradient
% simparams.optoptions = optimoptions('fmincon','Algorithm','interior-point','MaxFunctionEvaluations',3e5,'MaxIterations',1e4,'SpecifyConstraintGradient',true);
% simparams.optoptions = optimoptions('fmincon','Algorithm','sqp','MaxFunctionEvaluations',3e5,'MaxIterations',1e4,'SpecifyConstraintGradient',true);

%%%%%%%%%%%%%%%%%%% Objective and constraint function gradients 
simparams.optoptions = optimoptions('fmincon','Algorithm','interior-point','MaxFunctionEvaluations',3e5,'MaxIterations',1e4,'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);
% simparams.optoptions = optimoptions('fmincon','Algorithm','sqp','MaxFunctionEvaluations',3e5,'MaxIterations',1e4,'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);
% outputGradients = 1;
% outputCGradients = 1;

% Optimality and constraint satisfaction tolerances
simparams.optoptions.OptimalityTolerance = 1e-10;
simparams.optoptions.ConstraintTolerance = 1e-10 * simparams.m * simparams.n;
simparams.optoptions.StepTolerance = 1e-14; % use with sqp
% simparams.optoptions.FiniteDifferenceStepSize = 1e-5;

% To use parallel processing
% simparams.optoptions.UseParallel = true;

% Fmincon interior point feasibility mode
% simparams.optoptions.EnableFeasibilityMode = true;
% simparams.optoptions.SubproblemAlgorithm = 'cg';


% To have matlab check the analytical gradient, if being used
% simparams.optoptions = optimoptions('fmincon','Algorithm','interior-point','MaxFunctionEvaluations',3e5,'MaxIterations',1e4,'SpecifyObjectiveGradient',true,'CheckGradients',true,'FiniteDifferenceType','central','FiniteDifferenceStepSize',1e-8);


