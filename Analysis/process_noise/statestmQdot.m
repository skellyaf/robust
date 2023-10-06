function [xDot] = statestmQdot(t, X, mu, Qt)

% unpack x
assert(length(X) == 78);

G = [zeros(3,3); eye(3,3)];

x = X(1);
y = X(2);
z = X(3);
vx = X(4);
vy = X(5);
vz = X(6);


stm = reshape(X(7:42),6,6);
Q = reshape(X(43:78),6,6);

% X dot equations of motion for CR3BP
r1 = sqrt(  (x + mu)^2 + y^2 + z^2  );
r2 = sqrt(  (x-1+mu)^2 + y^2 + z^2  );

ax = 2*vy + x - (1-mu)*(x+mu)/r1^3 - mu*(x-1+mu)/r2^3;
ay = -2*vx + y -(1-mu)*y/r1^3 - mu*y/r2^3;
az = -(1-mu)*z/r1^3 - mu*z/r2^3;


F = cr3bp_sFrame_nd_Amatrix(X(1:6), mu);

stmdot = F*stm;

% Qdot equations of motion
Qdot = F*Q + Q*F' + G*Qt*G';




% Put it back together
xDot = [vx; vy; vz; ax; ay; az; stmdot(:); Qdot(:)];













end

