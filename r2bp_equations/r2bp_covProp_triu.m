function [xPdot] = r2bp_covProp(t, xP)
global mu;
tri_u = [1, 7, 8, 13, 14, 15, 19, 20, 21, 22, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36]';
diag_el = [1, 8, 15, 22, 29, 36];

% Parse inputs
x = xP(1:6);
P = zeros(6,6);
P(tri_u) = xP(7:end);
P = P' + P;
P(diag_el) = P(diag_el) / 2;
% P = reshape(xP(7:42),6,6);
r = [x(1); x(2); x(3)];
v = [x(4); x(5); x(6)];

% Two body dynamics
f = [zeros(3,3), eye(3,3); 
    -mu/norm(r)^3 * eye(3,3), zeros(3,3)];

xdot = f * x;


% Partial wrt state vector
A = r2bp_A_matrix(x);

% Cov prop, no noise
Pdot = A*P + P*A';

% Return
xPdot = [xdot; Pdot(tri_u)];
% xPdot = [xdot; reshape(Pdot,36,1)];

end