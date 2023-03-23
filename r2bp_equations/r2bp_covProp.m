function [xPdot] = r2bp_covProp(t, xP, mu)

% Parse inputs
x = xP(1:6);
P = reshape(xP(7:42),6,6);
r = [x(1); x(2); x(3)];
v = [x(4); x(5); x(6)];

% Two body dynamics
f = [zeros(3,3), eye(3,3); 
    -mu/norm(r)^3 * eye(3,3), zeros(3,3)];

xdot = f * x;


% Partial wrt state vector
A = r2bp_A_matrix(x, mu);

% Cov prop, no noise
Pdot = A*P + P*A';

% Return
xPdot = [xdot; reshape(Pdot,36,1)];

end