function [r_sun_sc, v_sc_sb1_nd] = calc_v_sc_lunarOrbit_br4bp(body, nu_moon, altitude, theta_em, simparams)

ndDist2km = simparams.l_star;;
ndVel2kms = simparams.l_star * simparams.n_sb1;

ndTime2sec = 1/simparams.n_sb1;

mub = simparams.mub;
mu = simparams.mu;

a4 = simparams.a4;


% Vectors to bodies
vecr4 = [-mub; 0; 0];

xe = 1 - mub - 1/a4 * mu * cos(theta_em);
ye = -1/a4 * mu *sin(theta_em);

xm = 1 - mub + 1/a4 * (1-mu) * cos(theta_em);
ym = 1/a4 * (1-mu) * sin(theta_em);

vecr1 = [xe; ye; 0];
vecr2 = [xm; ym; 0];






% r_earth = [xe; ye; 0];
r_moon = [xm; ym; 0];
r_body_sc_mag_km = altitude + body.rad;

r_moon_sc_km = [cos(nu_moon) * r_body_sc_mag_km; sin(nu_moon) * r_body_sc_mag_km; 0];
r_moon_sc = r_moon_sc_km / ndDist2km;
r_sun_sc = r_moon + r_moon_sc;

% Start with 2 body velocity relative to earth in km/s
v_sc_relMoon_kms = [-sin(nu_moon) * sqrt(body.mu/norm(r_moon_sc_km)); cos(nu_moon) * sqrt(body.mu/norm(r_moon_sc_km)); 0];
% Convert to ND
v_sc_relMoon_nd = v_sc_relMoon_kms / ndVel2kms;











% Calculate the velocity of Moon in the E-M frame (rotation around B1)
omega_em = [0; 0; simparams.theta_em_dot + 1];
% omega_em = [0; 0; n_em*ndTime2sec];
r_b1 = [1-mub; 0; 0];
r_b1_moon = r_moon - r_b1;
v_moon_relB1_nd = cross(omega_em, r_b1_moon);


% Calculate the velocity of B1 in the S-B1 frame
omega_sb1 = [0; 0; 1];
v_B1_relSun_nd = cross(omega_sb1, r_b1);


% Combine 
v_sc_sb1_nd = v_sc_relMoon_nd + v_moon_relB1_nd + v_B1_relSun_nd - cross(omega_sb1, r_sun_sc);





end