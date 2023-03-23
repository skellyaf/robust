%% Initial guess close to Hohmann for correction maneuver only


delta_t1 = 0;
x1 = [x_initial; delta_t1];

t_2_0 = t_0 + delta_t1;
% t_2_f = t_hohmann * .038; 
% t_2_f = t_hohmann * .0403195; 
% t_2_f = t_hohmann * .4; 
% t_2_f = t_hohmann * .17; 
% t_2_f = t_hohmann * .02; 
t_2_f = t_hohmann * .1; 
t_3_0 = t_2_f;

t_4_f = t_f;
t_4_0 = t_4_f - t_hohmann * .2;
t_3_f = t_4_0;

t_5_0 = t_4_f;
t_5_f = t_5_0;

delta_t2 = t_2_f-t_2_0;
% leg_t2 = delta_t2;
% leg_t2 = t_hohmann * .038;
% leg_t2 = 551.146440114601 + 0.0502706293529482;
x2 = [r_initial; v_xfer_1 * v_initial/norm(v_initial); delta_t2];





[~, x_3_0] = covProp(x2(1:6), P_initial, delta_t2,mu);
delta_t3 = t_3_f - t_3_0;
x3 = [x_3_0; delta_t3];

[~, x_4_0] = covProp(x3(1:6), P_initial, delta_t3,mu);
delta_t4 = t_4_f - t_4_0; 
x4 = [x_4_0; delta_t4];

delta_t5 = 0;
x5 = [x_target; delta_t5];

