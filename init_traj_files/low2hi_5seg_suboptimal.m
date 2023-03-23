t_back = T_initial * .25;


t_1 = T_initial * .15;
delta_t1 = t_1 - t_0;

x10 = stateStmProp(x_initial, -t_back, mu);

x1 = [x10; delta_t1 + t_back];

x1f = stateStmProp(x10, delta_t1 + t_back, mu);
x_initial = x10;
x0=x_initial;

r1f = x1f(1:3);
v1f = x1f(4:6);
uh = cross(r1f,v1f)/norm(cross(r1f,v1f));

delta_t2 = t_hohmann * .1;

r2f = [(alt_target+earth.rad)*.05; (alt_target+earth.rad)*.7; (alt_target+earth.rad)*.3];
uh = cross(r1f,r2f)/norm(cross(r1f,r2f));

[ v20, v2f ] = lambert (mu, delta_t2, r1f, r2f, uh);

x2 = [r1f; v20; delta_t2];

delta_t5 = t_hohmann * .3;
x5f = [x_target];
r5f = x5f(1:3);
v5f = x5f(4:6);


x50 = stateStmProp_fixed(x5f, -delta_t5, t_step, mu);

[a,e,i,W,w,nu] = rv2orbel(x50(1:3),x50(4:6),mu);
i = i+pi/6;
[r50, v50] = orbel2rv(a,e,i,W,w,nu,mu);

x5 = [r50; v50; delta_t5];

r4f = r50;
r40 = [(alt_target+earth.rad)*.05; (alt_target+earth.rad)*.3; (alt_target+earth.rad)*.5];



% r40 = r1f; %%%%
delta_t4 = t_hohmann * .4;
uh = cross(r40,r4f)/norm(cross(r40,r4f));


[v40, v4f] = lambert(mu, delta_t4, r40, r4f, uh)




uh = -cross(r50,r40)/norm(cross(r50,r40));

[ v40, v4f ] = lambert (mu, delta_t4, r40, r4f, uh);
x4 = [r40; v40; delta_t4];

r30 = r2f;

delta_t3 = t_hohmann - delta_t4 - delta_t5 - delta_t2 - delta_t1 + .5;

uh = -cross(r40,r30)/norm(cross(r40,r30));
[ v30, v3f ] = lambert (mu, delta_t3, r30, r40, uh);

x3 = [r30; v30; delta_t3];






t_corr = t_hohmann * .2;

x = [x1; x2; x3; x4; x5; t_corr]; 


x2h = [r_initial; v_xfer_1*v_initial/norm(v_initial)]
x3h = stateStmProp(x2h,t_hohmann/3, mu);
 
x4h = stateStmProp(x3h,t_hohmann/3, mu);

x_hohmann = [x10; t_back; 
    x2h; t_hohmann/3; 
    x3h; t_hohmann/3;
    x4h; t_hohmann/4;
    x_target; 0; 
    t_corr];

x_initial = x(1:6);
x0 = x_initial;
