function Xdot = br4bp_sb1_state_stm_stt_QdQ_de(t, X, mub, mu, a4, theta_em_dot, Qt)
%Function defining the differential equation for propagating the state
% state transition matrix, and state transition tensor in the BR4BP S-B1
% frame

G = [zeros(3,3); eye(3,3); zeros(1,3)];

nsv = 7;
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

F = br4bp_sb1_Amatrix(X(1:7), mub, mu, a4);
stmdot = F * stm;



stt2 = reshape(X(57:399), 7, 7, 7);


% STT propagation equation
Fab = br4bp_sb1_stt2_tensor(X(1:7), mub, mu, a4);
stt2dot = tensorCombine(stm,stt2,F,Fab);

% Q
Q = reshape(X(nsv+nsv^2+nsv^3+1:nsv+nsv^2+nsv^3+nsv^2),nsv,nsv);
% Qdot equations of motion
Qdot = F*Q + Q*F' + G*Qt*G';

% dQ
Q2 = reshape(X(nsv+2*nsv^2+nsv^3+1:nsv+2*nsv^2+2*nsv^3),nsv,nsv,nsv);

pFab = ipermute(Fab, [1, 2, 4, 3]); %T
pstm = ipermute(stm, [4, 3, 1, 2]); %A1

dFdx = sum(  permute(  pFab.*pstm, [1, 2, 3, 4]  ), [4]   );

Q2dot = tmult(F,Q2) +  tmult(Q2, F,[0 1]) + tmult(dFdx, Q) + tmult(Q, dFdx, [0 1]);




Xdot = [vx; vy; vz; ax; ay; az; theta_em_dot; stmdot(:); stt2dot(:); Qdot(:); Q2dot(:)];

end
