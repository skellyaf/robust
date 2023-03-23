function delr = cw_eqns(delx0,n,t)

% CW equations
% delx0 - state vector: [ delr0; delv0]
% n - mean motion
% t - time since initial position and velocity


delr0 = delx0(1:3);
delv0 = delx0(4:6);

delx = (4-3*cos(n*t)) * delr0(1) + sin(n*t)/n * delv0(1) + 2/n * (1-cos(n*t)) * delv0(2);

dely = 6*(sin(n*t)-n*t) * delr0(1) + delr0(2) + 2/n*(cos(n*t) - 1) * delv0(1) + 1/n*(4*sin(n*t)-3*n*t) * delv0(2);

delz = cos(n*t) * delr0(3) + 1/n*sin(n*t) * delv0(3);

delr = [delx; dely; delz];

end
