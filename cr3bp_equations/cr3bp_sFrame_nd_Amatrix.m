function A = cr3bp_sFrame_nd_Amatrix(X, mu)
%Function that returns the system matrix (A) given a state
% State - position and velocity of the negligible mass 3rd body in the
% circular restricted 3 body problem wrt the barycenter of the system (and
% frame that rotates with the earth-moon rotation rate)

% global muE muM Rm

assert(length(X)==6);

% Decompose states - position and velocity wrt barycenter of earth moon
R=X(1:3);
V=X(4:6);

x = R(1);
y = R(2);
z = R(3);

% mu = muM/(muE + muM);
r1 = sqrt(  (x + mu)^2 + y^2 + z^2  );
r2 = sqrt(  (x-1+mu)^2 + y^2 + z^2  );

% A matrix elements
% Build as a 6x6 matrix of 4 3x3 submatrices:
% RR RV
% VR VV

RR = zeros(3,3);
RV = eye(3);

dUdxx = 1 - (1-mu)/r1^3 - mu/r2^3 + 3*(1-mu)*(x+mu)^2/r1^5 + 3*mu*(x-1+mu)^2/r2^5;
dUdyy = 1 - (1-mu)/r1^3 - mu/r2^3 + 3*(1-mu)*y^2/r1^5 + 3*mu*y^2/r2^5;
dUdzz =   - (1-mu)/r1^3 - mu/r2^3 + 3*(1-mu)*z^2/r1^5 + 3*mu*z^2/r2^5;

dUdxy = 3*(1-mu)*(x+mu)*y/r1^5 + 3*mu*(x-1+mu)*y/r2^5;
dUdxz = 3*(1-mu)*(x+mu)*z/r1^5 + 3*mu*(x-1+mu)*z/r2^5;
dUdyz = 3*(1-mu)*y*z/r1^5 + 3*mu*y*z/r2^5;

VR = [dUdxx, dUdxy, dUdxz;
      dUdxy, dUdyy, dUdyz;
      dUdxz, dUdyz, dUdzz];

VV = zeros(3,3);
VV(1,2) = 2;
VV(2,1) = -2;

A = [RR RV; VR VV];


end

