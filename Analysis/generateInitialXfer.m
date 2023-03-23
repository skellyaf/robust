function [simparams] = generateInitialXfer(simparams)
%generateInitialXfer generates the initial n segments for a multi-segment
%transfer trajectory

mu = simparams.mu;

%% Generate an r and v vector from the orbital parameters in simparams
% Reshape from single trajectory vector into m x n segment matrix
simparams.x0 = reshape(simparams.x0,simparams.m, simparams.n);
% First segment
% Convert initial orbit COEs to R and V
coe = simparams.coe_init;

% Function accepts in meters, outputs in meters & m/s
if coe.inc == 0 && coe.ecc > 0 % equatorial elliptic, need to pass lonper
    [r0, v0] = keplerian2ijk(coe.a*1000, coe.ecc, coe.inc, coe.raan, coe.argp, coe.nu, 'longper', coe.lonper); % converting distance to m
elseif coe.inc == 0 && coe.ecc == 0 % equatorial circular, need to pass truelon
    [r0, v0] = keplerian2ijk(coe.a*1000, coe.ecc, coe.inc, coe.raan, coe.argp, coe.nu, 'truelon', coe.truelon); % converting distance to m
elseif coe.inc > 0 && coe.ecc == 0 % inclined circular, need to pass arglat
    [r0, v0] = keplerian2ijk(coe.a*1000, coe.ecc, coe.inc, coe.raan, coe.argp, coe.nu, 'arglat', coe.arglat); % converting distance to m
else % all the normal COEs apply
    [r0, v0] = keplerian2ijk(coe.a*1000, coe.ecc, coe.inc, coe.raan, coe.argp, coe.nu); % converting distance to m
end
simparams.x_init = [r0/1000; v0*3.6]; % currently converting back to km, km/hr
simparams.x0(1:6,1) = simparams.x_init;

% period of initial orbit
T0 = 2*pi*sqrt(simparams.coe_init.a^3/simparams.mu);
simparams.T0 = T0;
simparams.x0(7,1) = T0 * simparams.seg1_coast_fraction;

% Final segment
% Convert target orbit COEs to R and V
coe = simparams.coe_targ;

if coe.inc == 0 && coe.ecc > 0 % equatorial elliptic, need to pass lonper
    [rf, vf] = keplerian2ijk(coe.a*1000, coe.ecc, coe.inc, coe.raan, coe.argp, coe.nu, 'longper', coe.lonper); % converting distance to m
elseif coe.inc == 0 && coe.ecc == 0 % equatorial circular, need to pass truelon
    [rf, vf] = keplerian2ijk(coe.a*1000, coe.ecc, coe.inc, coe.raan, coe.argp, coe.nu, 'truelon', coe.truelon); % converting distance to m
elseif coe.inc > 0 && coe.ecc == 0 % inclined circular, need to pass arglat
    [rf, vf] = keplerian2ijk(coe.a*1000, coe.ecc, coe.inc, coe.raan, coe.argp, coe.nu, 'arglat', coe.arglat); % converting distance to m
else % all the normal COEs apply
    [rf, vf] = keplerian2ijk(coe.a*1000, coe.ecc, coe.inc, coe.raan, coe.argp, coe.nu); % converting distance to m
end
% [rf, vf] = keplerian2ijk(coe.a*1000, coe.ecc, coe.inc, coe.raan, coe.argp, coe.nu, 'truelon', coe.truelon);
% The target state:
simparams.x_target = [rf/1000; vf*3.6];

% period of target orbit
T_target = 2*pi*sqrt(simparams.coe_targ.a^3/simparams.mu);
simparams.T_target = T_target;
simparams.x0(7,simparams.n) = T_target * simparams.segn_coast_fraction;

% Propagate backwards to get state at the beginning of segment n
simparams.x0(1:6,simparams.n) = stateProp(simparams.x_target, -simparams.x0(end), simparams);

% Create n-2 intermediate segments
% First, use lambert to determine transfer initial state
% Duration/period of transfer - using period for a Hohmann xfer for now
a_xfer = (simparams.coe_init.a + simparams.coe_targ.a) / 2;
T_xfer = 2*pi*sqrt( (a_xfer)^3 / mu );
% DIVIDING HOHMANN XFER TIME BY 2 ---------
% T_xfer = 2*pi*1/sqrt( mu / norm(a_xfer)^3 )/2;

% Hohmann DV calc
v_xfer_1 = sqrt(2*mu/(norm(r0)/1000) - mu/a_xfer);
dv1 = v_xfer_1 - sqrt(mu/norm(r0/1000));
v_xfer_2 = sqrt(2*mu/(norm(rf)/1000) - mu/a_xfer);
dv2 = sqrt(mu/(norm(rf)/1000)) - v_xfer_2;
totalDVDet = dv1 + dv2

% Angular momentum
x20 = stateProp(simparams.x_init, simparams.x0(7), simparams);
r20 = x20(1:3);
r2f = simparams.x0(1:3,simparams.n);
uh = cross(r20,r2f)/norm(cross(r20,r2f));
if uh(3) < 0
    uh = -uh;
end

[v20, v2f] = lambert(mu, T_xfer/2, r20, r2f, uh);

% 2nd segment starts with results of lambert
simparams.x0(1:6,2) = [r20; v20];
simparams.x0(7,2) = T_xfer/2;

% Breaking the coast up into n-2 segments

if simparams.n >= 3
    % The number of intermediate coast segments
    n_int_seg = simparams.n-2;
    [~,x_t, t] = stateProp(simparams.x0(1:6,2), simparams.x0(7,2), simparams);

    % Dividing up the propagated states and times into n_int_seg segments
    idx_sep = ceil(length(t)/n_int_seg);

    % Segment indices from the propagation
    seg_idx = [1, (1:n_int_seg-1) * idx_sep];

    % Assign propagated coast states to segments
    simparams.x0(1:6,2:simparams.n-1) = x_t(seg_idx,1:6)';

    % Assign durations to each segment
    simparams.x0(7,2:end-2) = (t(seg_idx(2:end)) - t(seg_idx(1:end-1)))';
    simparams.x0(7,end-1) = t(end) - t(seg_idx(end));
    
end


% If want to use a fixed total transfer time, set it here:
simparams.tf = sum(simparams.x0(7,:));


% Reshape back to single vector
simparams.x0 = reshape(simparams.x0, simparams.n*simparams.m, 1);







%%
%% Saving how initialization used to happen below
%%


% 
% %% Define initial and target states
% % Target state in different altitude circular orbit
% alt_target = 35786;
% % alt_target = 1100;
% 
% % Initial state in circular orbit
% % 45 degrees inclined
% alt_0 = 450;
% r_initial = [alt_0 + earth.rad; 0; 0]; % All x position
% v_initial = [0; sqrt(mu/norm(r_initial)); 0]; % All y velocity
% 
% [a,e,i,W,w,nu] = rv2orbel(r_initial,v_initial,mu);
% i = i + pi/4;
% [r_initial, v_initial] = orbel2rv(a, e, i, 0, 0, 0, mu);
% 
% x_initial = [r_initial; v_initial];
% T_initial = 2*pi/sqrt( mu / norm(r_initial)^3 );
% 
% t_0 = 0; % initial time, zero HOURS
% 
% 
% 
% r_target = [- (alt_target + earth.rad); 0; 0]; % All -x position
% % inclined 45 degrees
% v_target = [0; -sqrt(mu/norm(r_target))*sqrt(2)/2; -sqrt(mu/norm(r_target))*sqrt(2)/2]; % All -y velocity 
% x_target = [r_target; v_target];
% T_target = 2*pi/sqrt( mu / norm(r_target)^3 );
% 
% 
% %% Define initial uncertainty
% 
% % load('C:\Users\skell\OneDrive - USU\Documents\code_repos\multiseg_opt\sims\20220708_1205.12_leo2leo_separate_EXAMPLE0p2\workspace.mat')
% 
% 
% 
% 
% %% Hohmann transfer orbit calculations - forms initial guess
% a_xfer = (norm(r_initial) + norm(r_target)) / 2;
% v_xfer_1 = sqrt(2*mu/norm(r_initial) - mu/a_xfer);
% v_xfer_2 = sqrt(2*mu/norm(r_target) - mu/a_xfer);
% 
% T_xfer = 2*pi*1/sqrt( mu / norm(a_xfer)^3 );
% t_hohmann = T_xfer/2;
% 
% 
% 
% 
% 
% %% Initial guess options
% 
% low2hi_3seg_suboptimal;
% % low2hi_3seg_suboptimal_2;
% % low2hi_5seg_suboptimal;
% % low2hi_closetohohmann;
% % low2hi_3seg_simple;
% 
% dv1dvctied = 0; % Flag to constrain the correction at the same time as DV1
% 
% x_hohmann = [x10; t_back; r_initial; v_xfer_1*v_initial/norm(v_initial); t_hohmann; x_target; 0; t_corr];
% 
% 
% 
% 
% % x = x_opt; % when loading a pre-computed trajectory
% %%
% 
% % A LEO to LEO poor initial guess
% % load('x_poor.mat');
% % x = x_poor(1:7,:);
% 
% 
% % removing the TCM time from the x vector
% 
% x = x(1:21);
% x = x_hohmann(1:21);

end