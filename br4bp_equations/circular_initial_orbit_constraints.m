function [ceq, ceq_grad] = circular_initial_orbit_constraints(x0, r_e, fixed_orbit_radius_nd, mu_earth_nd, simparams)
%circular_initial_orbit_constraints Summary of this function goes here
%   Detailed explanation goes here

ceq = zeros(3,1);
ceq_grad = zeros(3,8);

% Setup
x0_1 = x0;
theta_em01 = x0(7);
mub = simparams.mub;
mu = simparams.mu;
a4 = simparams.a4;


% Passing these in instead of calculating
% xe_0 = 1 - mub - 1/a4 * mu * cos(theta_em01);
% ye_0 = -1/a4 * mu *sin(theta_em01);
% r_e = [xe_0; ye_0; 0];




%%
%%%%%%%% Constraint 1 - position wrt body magnitude (distance from body)
r_sc = x0_1(1:3);
r_sc_e = r_sc - r_e;
r_sc_e_mag = norm(r_sc_e);

% Constraint:
% r_fixed_e = body.rad + fixed_orbit_altitude_nd;
ceq(1) = r_sc_e_mag - fixed_orbit_radius_nd; % in ND distance

% Gradient
% drscE_dr01
i_rscE = r_sc_e / r_sc_e_mag;
drscEmag_dr01 = i_rscE';
% drscE_dv01
drscEmag_dv01 = zeros(1,3);
% drscE_dth01
drE_dth01 = [mu/a4 * sin(theta_em01); -mu/a4 * cos(theta_em01); 0];
drscEmag_dth01 = - i_rscE' * drE_dth01;
% drscE_ddt1
drscEmag_ddt1 = 0;
% Putting it together
ceq_grad(1,:) = [drscEmag_dr01, drscEmag_dv01, drscEmag_dth01, drscEmag_ddt1];



%%
%%%%%%% Constraint 2 - fixed 2 body orbital energy wrt earth
% [r_sun_sc, v_sc_sb1] = calc_v_sc_earthOrbit_br4bp(body, nu_earth, altitude, theta_em01, simparams);
v_sc_sb1 = x0_1(4:6);

% mu_earth_nd = earth.mu / ndDist2km^3 * ndTime2sec^2;

% v_EB1_B1
w_em = [0; 0; simparams.theta_em_dot + 1];
w_em_x = cross_prod(w_em);

r_B1 = [1-mub; 0; 0];
r_EB1 = r_e - r_B1;
v_EB1_B1 = w_em_x * r_EB1;

% v_B1S_S
w_sB1 = [0; 0; 1];
w_sB1_x = cross_prod(w_sB1);
r_B1S = r_B1;
v_B1S_S = w_sB1_x * r_B1S;

r_scS = r_sc;
v_scE_E = v_sc_sb1 - v_EB1_B1 - v_B1S_S + w_sB1_x * r_scS;
v_scE_E_mag = norm(v_scE_E);
eps_0 = v_scE_E_mag^2 / 2 - mu_earth_nd / r_sc_e_mag;
% Verification
% r = earth.rad + altitude;
r = fixed_orbit_radius_nd;
% v = sqrt(mu_earth_nd/r);
% eps_fixed_km2s2 = v^2/2 - mu_earth_nd/r;
% eps_fixed_0 = eps_fixed_km2s2 / ndVel2kms^2;
eps_fixed_0 = -mu_earth_nd / (2*r);

% Constraint 2:
ceq(2) = eps_0 - eps_fixed_0;

% Gradient
i_v_scE_E = v_scE_E / v_scE_E_mag;

% dvscSB1_ds1
dvscSB1_ds1 = [zeros(3,3), eye(3,3), zeros(3,2)];

% dvEB1B1_ds1
drEB1_dth01 = drE_dth01;
dvEB1B1_ds1 = [zeros(3,6), w_em_x * drEB1_dth01, zeros(3,1)];

% dvB1SS_ds1
dvB1SS_ds1 = zeros(3,8);

% dwSB1xrsc_ds1
dwSB1xrsc_ds1 = w_sB1_x * [eye(3), zeros(3,5)];

% Sum them
dvscEE_ds1 = dvscSB1_ds1 - dvEB1B1_ds1 - dvB1SS_ds1 + dwSB1xrsc_ds1;
dvscEEmag_ds1 = i_v_scE_E' * dvscEE_ds1;


% They contribute to the gradient of the energy equation
deps0_dr01 = mu_earth_nd / r_sc_e_mag^2 * i_rscE';
deps0_dv01 = v_scE_E' * dvscEE_ds1(:,4:6);
deps0_dth01 = v_scE_E' * dvscEE_ds1(:,7) + mu_earth_nd / r_sc_e_mag^2 * drscEmag_dth01;
deps0_ddt1 = 0;
deps0_ds1 = [deps0_dr01, deps0_dv01, deps0_dth01, deps0_ddt1];
ceq_grad(2,:) = deps0_ds1;



%%
%%%%%%% Constraint 3 - orthogonal r and v vectors
ceq(3) = r_sc_e' * v_scE_E;


% Gradient
drscE_dr01 = eye(3);
drscE_dv01 = zeros(3,3);
drscE_dth01 = - drE_dth01;
drscE_ddt1 = zeros(3,1);

drscE_ds1 = [drscE_dr01, drscE_dv01, drscE_dth01, drscE_ddt1];

% dvscEE_ds1 already exists

drtvE_ds1 = v_scE_E' * drscE_ds1 + r_sc_e' * dvscEE_ds1;
ceq_grad(3,:) = drtvE_ds1;


end