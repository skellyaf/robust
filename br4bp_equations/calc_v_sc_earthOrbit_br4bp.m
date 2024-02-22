function [r_sun_sc, v_sc_sb1_nd] = calc_v_sc_earthOrbit_br4bp(body, nu_earth, altitude, theta_em0, simparams)

ndDist2km = simparams.l_star;
ndVel2kms = simparams.l_star * simparams.n_sb1;

mub = simparams.mub;
mu = simparams.mu;

a4 = simparams.a4;


% Vectors to bodies
vecr4 = [-mub; 0; 0];

xe = 1 - mub - 1/a4 * mu * cos(theta_em0);
ye = -1/a4 * mu *sin(theta_em0);

xm = 1 - mub + 1/a4 * (1-mu) * cos(theta_em0);
ym = 1/a4 * (1-mu) * sin(theta_em0);

vecr1 = [xe; ye; 0];
vecr2 = [xm; ym; 0];






r_earth = [xe; ye; 0];
r_earth_sc_mag_km = altitude + body.rad;

r_earth_sc_km = [cos(nu_earth) * r_earth_sc_mag_km; sin(nu_earth) * r_earth_sc_mag_km; 0];
r_earth_sc = r_earth_sc_km / ndDist2km;
r_sun_sc = r_earth + r_earth_sc;

% Start with 2 body velocity relative to earth in km/s
v_sc_relEarth_kms = [-sin(nu_earth) * sqrt(body.mu/norm(r_earth_sc_km)); cos(nu_earth) * sqrt(body.mu/norm(r_earth_sc_km)); 0];
% Convert to ND
v_sc_relEarth_nd = v_sc_relEarth_kms / ndVel2kms;











% Calculate the velocity of Earth in the E-M frame (rotation around B1)
omega_em = [0; 0; simparams.theta_em_dot + 1];
% omega_em = [0; 0; n_em*ndTime2sec];
r_b1 = [1-mub; 0; 0];
r_b1_earth = r_earth - r_b1;
v_earth_relB1_nd = cross(omega_em, r_b1_earth);


% Calculate the velocity of B1 in the S-B1 frame
omega_sb1 = [0; 0; 1];
v_B1_relSun_nd = cross(omega_sb1, r_b1);


% Combine 
v_sc_sb1_nd = v_sc_relEarth_nd + v_earth_relB1_nd + v_B1_relSun_nd - cross(omega_sb1, r_sun_sc);





end