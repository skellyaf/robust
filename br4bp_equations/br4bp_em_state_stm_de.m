function Xdot = br4bp_em_state_stm_de(t, X, mu, m4, a4, theta_s_dot)
%br4bp_em_state_stm_de is the first-order time differential equations of
%motion for the bicircular restricted 4 body problem in the earth-moon
%rotating frame.
% X is the 7 dimensional state vector [r, v, theta_s]
% theta_s is the sun angle

% Unpack state
x = X(1);
y = X(2);
z = X(3);
vx = X(4);
vy = X(5);
vz = X(6);
theta_s = X(7);


r1 = sqrt(  (x + mu)^2 + y^2 + z^2  );
r2 = sqrt(  (x-1+mu)^2 + y^2 + z^2  );
r4 = sqrt( (x - a4 * cos(theta_s))^2 + (y - a4 * sin(theta_s))^2 + z^2 );


ax = 2*vy + x - (1-mu)*(x+mu)/r1^3 - mu*(x-1+mu)/r2^3 ...
     - m4 / r4^3 * (x - a4 * cos(theta_s)) - m4 * cos(theta_s) / a4^2;

ay = -2*vx + y -(1-mu)*y/r1^3 - mu*y/r2^3 ...
     - m4 / r4^3 * (y - a4 * sin(theta_s)) - m4 * sin(theta_s) / a4^2;

az = -(1-mu)*z/r1^3 - mu*z/r2^3 ...
     - m4 / r4^3 * z;



stm = reshape(X(8:56), 7, 7);

A = br4bp_em_Amatrix(X(1:7), mu, m4, a4);
stmdot = A * stm;



Xdot = [vx; vy; vz; ax; ay; az; theta_s_dot; reshape(stmdot, 49, 1)];




end