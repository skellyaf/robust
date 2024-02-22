function Xdot = br4bp_sb1_state_stm_stt_de(t, X, mub, mu, a4, theta_em_dot)
%Function defining the differential equation for propagating the state
% state transition matrix, and state transition tensor in the BR4BP S-B1
% frame

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


stm = reshape(X(8:56), 7, 7);

A = br4bp_sb1_Amatrix(X(1:7), mub, mu, a4);
stmdot = A * stm;



stt2 = reshape(X(57:399), 7, 7, 7);


% STT propagation equation
Fab = br4bp_sb1_stt2_tensor(X(1:7), mub, mu, a4);
stt2dot = tensorCombine(stm,stt2,A,Fab);

Xdot = [vx; vy; vz; ax; ay; az; theta_em_dot; stmdot(:); stt2dot(:)];

end
