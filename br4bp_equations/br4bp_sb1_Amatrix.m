function A = br4bp_sb1_Amatrix(X, mub, mu, a4)
%br4bp_sb1_Amatrix returns the system matrix in the BR4BP Sun-B1 frame
% Inputs are the mass ratio of the Sun-Earth-Moon system (mub =
% mu_underbar), the mass ratio of the Earth-Moon system (mu), and the
% semimajor axis of the Earth's orbit around the sun (a4).

assert(length(X)==7);

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

% r1 = sqrt(  (x + mu)^2 + y^2 + z^2  );
% r2 = sqrt(  (x-1+mu)^2 + y^2 + z^2  );

% r4 = sqrt( (x - a4 * cos(theta_s))^2 + (y - a4 * sin(theta_s))^2 + z^2 );

% A matrix elements

A = zeros(7,7);
A(1:3,4:6) = eye(3);

dYdxx = 1 - (1-mub)/r43^3 + 3*(1-mub)*(x+mub)^2/r43^5 - mub*(1-mu)/r13^3 + 3*mub*(1-mu)*(x-xe)^2/r13^5 + ...
        - mub*mu/r23^3 + 3*mub*mu*(x-xm)^2/r23^5; 

dYdyy = 1 - (1-mub)/r43^3 + 3*(1-mub)*y^2/r43^5 - mub*(1-mu)/r13^3 + 3*mub*(1-mu)*(y-ye)^2/r13^5 + ...
        - mub*mu/r23^3 + 3*mub*mu*(y-ym)^2/r23^5; 

dYdzz = - (1-mub)/r43^3 + 3*(1-mub)*z^2/r43^5 - mub*(1-mu)/r13^3 + 3*mub*(1-mu)*z^2/r13^5 + ...
         - mub*mu/r23^3 + 3*mub*mu * z^2 / r23^5;

dYdxy = 3*(1-mub)*(x+mub)*y / r43^5 + 3*mub*(1-mu)*(x-xe)*(y-ye) / r13^5 + ...
        3*mub*mu*(x-xm)*(y-ym) / r23^5;


dYdxz = 3*(1-mub)*(x+mub)*z / r43^5 + 3*mub*(1-mu)*(x-xe)*z / r13^5 + ...
        3*mub*mu*(x-xm)*z / r23^5;


dYdyz = 3*(1-mub)*y*z / r43^5 + 3*mub*(1-mu)*(y-ye)*z / r13^5 + ...
        3*mub*mu*(y-ym)*z / r23^5;




cse1 = (y-ye)*cos(theta_em) - (x-xe)*sin(theta_em);
scm2 = (x-xm)*sin(theta_em) - (y-ym)*cos(theta_em);

dYdxdtheta = 3*mub*mu*(1-mu)*(x-xe)*cse1 / r13^5 / a4 ...
           + mub*mu*(1-mu)*sin(theta_em) / r13^3 / a4 ...
           + 3*mub*mu*(1-mu)*(x-xm)*scm2 / r23^5 / a4 ...
           - mub*mu*(1-mu)*sin(theta_em) / r23^3 / a4;

dYdydtheta = 3*mub*mu*(1-mu)*(y-ye)*cse1 / r13^5 / a4 ...
           - mub*mu*(1-mu)*cos(theta_em) / r13^3 / a4 ...
           + 3*mub*mu*(1-mu)*(y-ym)*scm2 / r23^5 / a4 ...
           + mub*mu*(1-mu)*cos(theta_em) / r23^3 / a4;

dYdzdtheta = 3*mub*mu*(1-mu)*(z)*cse1 / r13^5 / a4 ...
           + 3*mub*mu*(1-mu)*(z)*scm2 / r23^5 / a4;




A(4:6,1:3) = [dYdxx, dYdxy, dYdxz;
              dYdxy, dYdyy, dYdyz;
              dYdxz, dYdyz, dYdzz];

A(4:6,7) = [dYdxdtheta; dYdydtheta; dYdzdtheta];

A(4,5) = 2;
A(5,4) = -2;



end