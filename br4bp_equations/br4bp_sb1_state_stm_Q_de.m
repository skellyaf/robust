function Xdot = br4bp_sb1_state_stm_Q_de(t, X, mub, mu, a4, theta_em_dot, Qt)
%br4bp_em_state_stm_de is the first-order time differential equations of
%motion for the bicircular restricted 4 body problem in the sun-B1
%rotating frame.
% X is the 7 dimensional state vector [r, v, theta_em]
% theta_em is the Earth-Moon line angle measured from B1 along the x-axis
% positive clockwise.

nsv = 7;
assert(length(X) == 7 + 7*7 + 7*7);
G = [zeros(3,3); eye(3,3); zeros(1,3)];

% Unpack state
x = X(1);
y = X(2);
z = X(3);
vx = X(4);
vy = X(5);
vz = X(6);
theta_em = X(7);

vecr3 = [x; y; z];
vecr4 = [-mub; 0; 0];

xe = 1 - mub - 1/a4 * mu * cos(theta_em);
ye = -1/a4 * mu *sin(theta_em);

xm = 1 - mub + 1/a4 * (1-mu) * cos(theta_em);
ym = 1/a4 * (1-mu) * sin(theta_em);

vecr1 = [xe; ye; 0];
vecr2 = [xm; ym; 0];

vecr13 = vecr3 - vecr1;
vecr23 = vecr3 - vecr2;
vecr43 = vecr3 - vecr4;

r13 = vecnorm(vecr13);
r23 = vecnorm(vecr23);
r43 = vecnorm(vecr43);


dYdx = x - (1-mub)*(x+mub) / r43^3 - mub*(1-mu)*(x-xe) / r13^3 ...
     - mub*mu*(x-xm) / r23^3;

dYdy = y - (1-mub)*y / r43^3 - mub*(1-mu)*(y-ye) / r13^3 ...
     - mub*mu*(y-ym) / r23^3;

dYdz = - (1-mub)*z / r43^3 - mub*(1-mu)*z / r13^3 ...
     - mub*mu*z / r23^3;


ax = 2*vy + dYdx;
ay = -2*vx + dYdy;
az = dYdz;


stm = reshape(X(nsv + 1:nsv + nsv*nsv), nsv, nsv);

F = br4bp_sb1_Amatrix(X(1:nsv), mub, mu, a4);
stmdot = F * stm;

Q = reshape(X(nsv + nsv*nsv + 1:nsv + nsv*nsv + nsv*nsv),nsv,nsv);


% Qdot equations of motion
Qdot = F*Q + Q*F' + G*Qt*G';

% Put it back together
Xdot = [vx; vy; vz; ax; ay; az; theta_em_dot; stmdot(:); Qdot(:)];




end