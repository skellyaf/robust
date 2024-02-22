function A = br4bp_em_Amatrix(X, mu, m4, a4)
%br4bp_em_Amatrix returns the system matrix in the BR4BP given a state and
%mass of the sun

assert(length(X)==7);

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

% A matrix elements

A = zeros(7,7);
A(1:3,4:6) = eye(3);

dYdxx = 1 - (1-mu)/r1^3 - mu/r2^3 + 3*(1-mu)*(x+mu)^2/r1^5 + 3*mu*(x-1+mu)^2/r2^5 + ...
        3 * m4/r4^5 * (x - a4*cos(theta_s))^2 - m4/r4^3;

dYdyy = 1 - (1-mu)/r1^3 - mu/r2^3 + 3*(1-mu)*y^2/r1^5 + 3*mu*y^2/r2^5 + ...
        3 * m4/r4^5 * (y - a4*sin(theta_s))^2 - m4/r4^3;

dYdzz = - (1-mu)/r1^3 - mu/r2^3 + 3*(1-mu)*z^2/r1^5 + 3*mu*z^2/r2^5 + ...
        3 * m4 / r4^5 * z^2 - m4/r4^3;

dYdxy = 3*(1-mu)*(x+mu)*y/r1^5 + 3*mu*(x-1+mu)*y/r2^5 + ...
        3 * m4/r4^5 * (x - a4*cos(theta_s)) * (y - a4*sin(theta_s));
dYdxz = 3*(1-mu)*(x+mu)*z/r1^5 + 3*mu*(x-1+mu)*z/r2^5 + ...
        3 * m4/r4^5 * z * (x - a4*cos(theta_s));
dYdyz = 3*(1-mu)*y*z/r1^5 + 3*mu*y*z/r2^5 + ...
        3 * m4/r4^5 * z * (y - a4*sin(theta_s));

dYdxdtheta = 3/2 * m4/r4^5 * (x-a4*cos(theta_s)) * (2*x*a4*sin(theta_s) - a4^2*sin(2*theta_s)) - ...
             m4/r4^3 * a4*sin(theta_s) + m4/a4^2 * sin(theta_s);

dYdydtheta = 3/2 * m4/r4^5 * (x-a4*sin(theta_s)) * (2*x*a4*sin(theta_s) - a4^2*sin(2*theta_s)) + ...
             m4/r4^3 * a4*cos(theta_s) - m4/a4^2 * cos(theta_s);

dYdzdtheta = 3/2 * m4/r4^5 * z * (2*x*a4*sin(theta_s) - a4^2*sin(2*theta_s));

A(4:6,1:3) = [dYdxx, dYdxy, dYdxz;
      dYdxy, dYdyy, dYdyz;
      dYdxz, dYdyz, dYdzz];

A(4:6,7) = [dYdxdtheta; dYdydtheta; dYdzdtheta];

A(4,5) = 2;
A(5,4) = -2;



end