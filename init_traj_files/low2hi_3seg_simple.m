% Hohmann transfer orbit - initial guess
a_xfer = (norm(r_initial) + norm(r_target)) / 2;
v_xfer_1 = sqrt(2*mu/norm(r_initial) - mu/a_xfer);

T_xfer = 2*pi*1/sqrt( mu / norm(a_xfer)^3 );
t2_hohmann = T_xfer/2;

%% Propagate initial state backward
% To simulate a coast into the initial orbit prior to Hohmann xfer


% Segment 1 - propagate backwards from starting position
t1_hohmann = T_initial/4; % 10 minute coast 
% Propagating initial orbit backwards to get starting state
[t, x_initial_back] = ode45(@r2bp_de, [0,-t1_hohmann], x_initial, mu, options);
x0 = x_initial_back(end,:)';

t_f = t1_hohmann + t2_hohmann;
t1 = t1_hohmann;
delta_t1=t1_hohmann;

%% Propagate target state backwards also

delta_t3 = T_target * .125;
[t, x_3_back] = ode45(@r2bp_de, [0,-delta_t3], x_target, mu, options);
x_3_0 = x_3_back(end,:)';


delta_t2 = t2_hohmann - delta_t3;

x3 = [x_3_0; delta_t3];



%% Segment definition
% Segment 1 state definition
% x1 = [x0; reshape(P_initial,36,1); t_0; t_0 + t1];
x1 = [x0; delta_t1];

% Segment 2
v_1 = [0; v_xfer_1; 0];
% Propagate covariance over initial coast ( not required in rev 2
% [P_1_end, x_1_end, x_1_t] = covProp(x1(1:6), P_initial, delta_t1);

% x2 = [r_initial; v_1; reshape(P_1_end,36,1); t_0 + delta_t1; t_0 + delta_t1 + delta_t2];


% Lambert targeting to get from r20 to r30
r20 = r_initial;
r30 = x_3_0(1:3);
uh = cross(r20,r30)/norm(cross(r20,r30));
[v20, v2f] = lambert(mu,delta_t2,r20, r30,uh);



x2 = [r20; v20; delta_t2];
% x2 = [r_initial; v_1; t_0 + delta_t1; t_0 + delta_t1 + delta_t2];

% Segment 3
% Propagate covariance over transfer orbit - not required in rev 2
% [P_2_end, x_2_end, x_2_t] = covProp(x2(1:6), P_1_end, delta_t2);
% Add delta V to transfer to circularize at higher orbit
% deltaV2 = norm(v_target) - norm(x_2_end(4:6));
% v3 = x_2_end(4:6) + deltaV2 * x_2_end(4:6)/norm(x_2_end(4:6));
% x3 = [x_2_end(1:3); v3; reshape(P_2_end,36,1); t_0 + delta_t1 + delta_t2; t_0 + delta_t1 + delta_t2 + delta_t3];
% x3 = [x_2_end(1:3); v3; t_0 + delta_t1 + delta_t2; t_0 + delta_t1 + delta_t2 + delta_t3];


% End of Segment 3 to observe initial guess covariance at end of sim
% [P_3_end, x_3_end, x_3_t] = covProp(x3(1:6), P_2_end, delta_t3);
% P_3_end = P_2_end;
% trace(P_3_end(1:3,1:3))



% Verifying that STM prop is equivalent covariance to previous cov prop eqns
% Propagate STM for segment 1 from x_initial for t1_hohmann seconds
% [xfinal1, stm1] = stateStmProp(x0, delta_t1)
% [xfinal2, stm2] = stateStmProp(x2(1:6), delta_t2)
% [xfinal3, stm3] = stateStmProp(x3(1:6), delta_t3)
% Pfinal = stm2 * stm1 * P_initial * (stm2 * stm1)'
% trace(Pfinal(1:3,1:3))


% Issue with extra time - covariance in target orbit will grow also
% Think about how/when to include covariance measurement updates

x = [x1; x2; x3; 0];
