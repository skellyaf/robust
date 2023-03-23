%% Hohmann transfer orbit calculations - forms initial guess
a_xfer = (norm(r_initial) + norm(r_target)) / 2;
v_xfer_1 = sqrt(2*mu/norm(r_initial) - mu/a_xfer);
v_xfer_2 = sqrt(2*mu/norm(r_target) - mu/a_xfer);

T_xfer = 2*pi*1/sqrt( mu / norm(a_xfer)^3 );
t_hohmann = T_xfer/2;



%% 3 segment poor initial guess


t_1 = T_initial * .15;
delta_t1 = t_1 - t_0;
x1 = [x_initial; delta_t1];

x1f = stateStmProp(x_initial, delta_t1, mu);

r1f = x1f(1:3);
v1f = x1f(4:6);

delta_t2 = t_hohmann * .9;



delta_t3 = t_hohmann * .3;
x3f = [x_target];
r3f = x3f(1:3);
v3f = x3f(4:6);


x30 = stateStmProp(x3f, -delta_t3, mu);

[a,e,i,W,w,nu] = rv2orbel(x30(1:3),x30(4:6),mu);
i = i+pi/6;
[r30, v30] = orbel2rv(a,e,i,W,w,nu,mu);

x3 = [r30; v30; delta_t3];

r2f = r30;

r20 = r1f;
uh = cross(r20,r2f)/norm(cross(r20,r2f));


[v20, v2f] = lambert(mu, delta_t2, r20, r2f, uh);


x2 = [r20; v20; delta_t2];
x3 = [r30; v30; delta_t3];


% Propagate x1 backwards
T1 = 2*pi/sqrt(mu/norm(x1(1:3))^3);
% t_back = .2;
t_back = .25*T1;
t_f = t_hohmann + t_back;

x1 = stateStmProp(x1(1:6), -t_back, mu)
x1 = [x1; t_back + delta_t1];
x10 = x1(1:6);

t_corr = t_hohmann * .2 + t_back;
x = [x1; x2; x3; t_corr];
x_initial = x(1:6);
x0 = x_initial;
