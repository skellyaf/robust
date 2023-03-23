function Xdot = cr3bp_sFrame_nd_stt2_de(t, X, mu)
%Function defining the differential equation for propagating the state
% and state transition matrix in the circular restricted 3 body problem
% STM differential equation: stmdot = A*stm_last

% global muE muM

x = X(1);
y = X(2);
z = X(3);
vx = X(4);
vy = X(5);
vz = X(6);

stm = reshape(X(7:42), 6,6);

% ADDED STT
stt2 = reshape(X(43:258),6,6,6);

% mu = muM/(muE + muM);

r1 = sqrt(  (x + mu)^2 + y^2 + z^2  );
r2 = sqrt(  (x-1+mu)^2 + y^2 + z^2  );

ax = 2*vy + x - (1-mu)*(x+mu)/r1^3 - mu*(x-1+mu)/r2^3;
ay = -2*vx + y -(1-mu)*y/r1^3 - mu*y/r2^3;
az = -(1-mu)*z/r1^3 - mu*z/r2^3;




A = cr3bp_sFrame_nd_Amatrix(X(1:6), mu);
stmdot = A*stm;



% STT propagation equation
Fab = cr3bp_stt2_tensor([x; y; z; vx; vy; vz], mu);
stt2dot = tensorCombine(stm,stt2,A,Fab);

Xdot = [vx; vy; vz; ax; ay; az; stmdot(:); stt2dot(:)];


end
