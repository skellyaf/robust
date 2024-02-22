function [ceq, ceq_grad] = circular_target_orbit_constraints(xf_n, r_mf, stm_n, fixed_orbit_radius_nd, mu_moon_nd, simparams)
%circular_target_orbit_constraints Summary of this function goes here
%   Detailed explanation goes here

ceq = zeros(3,1);
ceq_grad = zeros(3,8);

% Setup
theta_emfn = xf_n(7);
mub = simparams.mub;
mu = simparams.mu;
a4 = simparams.a4;

% Passing in the body position instead of calculating it inside
% xm_f = 1 - mub + 1/a4 * (1-mu) * cos(theta_emfn);
% ym_f = 1/a4 * (1-mu) * sin(theta_emfn);
% 
% r_mf = [xm_f; ym_f; 0];

r_B1 = [1-mub; 0; 0];
w_em = [0; 0; simparams.theta_em_dot + 1];
w_em_x = cross_prod(w_em);


%%%%%% Constraint 4 - C1 for the moon
r_sc = xf_n(1:3);
r_sc_m = r_sc - r_mf;
r_sc_m_mag = norm(r_sc_m);


% lunar_alt = 100;
% r_fixed_m = (moon.rad + lunar_alt) / ndDist2km;

ceq(1) = r_sc_m_mag - fixed_orbit_radius_nd;


% Gradient
% drscMmag_dr0n
i_rscM = r_sc_m / r_sc_m_mag;
drscMmag_dr0n = i_rscM' * stm_n(1:3,1:3);

% drscMmag_dV0n
drscMmag_dV0n = i_rscM' * stm_n(1:3,4:6);

% drscMmag_dth0n
drfn_dth0n = stm_n(1:3,7);

drMfn_dth0n = [-1/a4*(1-mu) * sin(theta_emfn); ...
               1/a4* (1-mu) * cos(theta_emfn); ...
               0];

drscMmag_dth0n = i_rscM' * (drfn_dth0n - drMfn_dth0n);

% drscMmag_ddtn
vf_n = xf_n(4:6);
drfn_ddtn = vf_n;
drscMmag_ddtn = i_rscM' * (drfn_ddtn - simparams.theta_em_dot * drMfn_dth0n);

drscMmag_dsn = [drscMmag_dr0n, drscMmag_dV0n, drscMmag_dth0n, drscMmag_ddtn];

ceq_grad(1,:) = drscMmag_dsn;



%%
%%%%%%% Constraint 5 - fixed 2 body orbital energy wrt moon 
% lunar_altitude = 100; % km, arrival altitude
% nu_moon = 270*pi/180;

% [r_sun_sc_f, v_sc_sb1_f] = calc_v_sc_lunarOrbit_br4bp(moon, nu_moon, lunar_altitude, theta_emfn, simparams);
v_sc_sb1_f = xf_n(4:6);
% mu_moon_nd = moon.mu / ndDist2km^3 * ndTime2sec^2;



% v_MB1_B1
% w_em = [0; 0; simparams.theta_em_dot + 1];
% w_em_x = cross_prod(w_em);

% r_B1 = [1-mub; 0; 0];
r_MB1 = r_mf - r_B1;
v_MB1_B1 = w_em_x * r_MB1;

% v_B1S_S
w_sB1 = [0; 0; 1];
w_sB1_x = cross_prod(w_sB1);
r_B1S = r_B1;
v_B1S_S = w_sB1_x * r_B1S;

r_scS = r_sc;
v_scM_M = v_sc_sb1_f - v_MB1_B1 - v_B1S_S + w_sB1_x * r_scS;
v_scM_M_mag = norm(v_scM_M);
eps_f = v_scM_M_mag^2 / 2 - mu_moon_nd / r_sc_m_mag;
% Verification
% r = moon.rad + lunar_altitude;
r = fixed_orbit_radius_nd;
% v = sqrt(moon.mu/r);
% eps_fixed_km2s2 = v^2/2 - moon.mu/r;
% eps_fixed_f = eps_fixed_km2s2 / ndVel2kms^2;
eps_fixed_f = -mu_moon_nd / (2 * r);

% Constraint 2:
ceq(2) = eps_f - eps_fixed_f;


% Gradient
i_v_scM_M = v_scM_M / v_scM_M_mag;

% dvscSB1_dsn
dvscSB1_dsn = [stm_n(4:6,:), vf_n];

% dvMB1B1_dsn
drMB1_dthfn = drMfn_dth0n;
dvMB1B1_dsn = w_em_x * [zeros(3,6), drMB1_dthfn, simparams.theta_em_dot * drMfn_dth0n];

% dvB1SS_dsn
dvB1SS_dsn = zeros(3,8);

% dwSB1xrsc_dsn
dwSB1xrsc_dsn = w_sB1_x * [stm_n(1:3,:), vf_n];

% Sum them
dvscMM_dsn = dvscSB1_dsn - dvMB1B1_dsn - dvB1SS_dsn + dwSB1xrsc_dsn;
dvscMMmag_dsn = i_v_scM_M' * dvscMM_dsn;


% They contribute to the gradient of the energy equation
depsf_dr0n = mu_moon_nd / r_sc_m_mag^2 * drscMmag_dr0n;
depsf_dv0n = v_scM_M' * dvscMM_dsn(:,4:6);
depsf_dth0n = v_scM_M' * dvscMM_dsn(:,7) + mu_moon_nd / r_sc_m_mag^2 * drscMmag_dth0n;
depsf_ddtn = v_scM_M' * dvscMM_dsn(:,8) + mu_moon_nd / r_sc_m_mag^2 * drscMmag_ddtn;
depsf_dsn = [depsf_dr0n, depsf_dv0n, depsf_dth0n, depsf_ddtn];

ceq_grad(2,:) = depsf_dsn;

%%
%%%%%%% Constraint 6 - orthogonal r and v vectors

ceq(3) = r_sc_m' * v_scM_M;


% Gradient
drscM_dr0n = stm_n(1:3,1:3);
drscM_dv0n = stm_n(1:3,4:6);
drscM_dth0n = stm_n(1:3,7) - drMfn_dth0n;
drscM_ddtn = drfn_ddtn - simparams.theta_em_dot * drMfn_dth0n;

drscM_dsn = [drscM_dr0n, drscM_dv0n, drscM_dth0n, drscM_ddtn];

% dvscMM_dsn already exists

drtvM_dsn = v_scM_M' * drscM_dsn + r_sc_m' * dvscMM_dsn;

ceq_grad(3,:) = drtvM_dsn;

end