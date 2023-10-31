clear; close all; clc;


%% Traj 1 - 12 hour trail
% C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\20231004_1516.52_rdvz_nrho_2dv_12hrlead_det_Qem8_4TcmMax_100mPmaxr

load('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\20231004_1516.52_rdvz_nrho_2dv_12hrlead_det_Qem8_4TcmMax_100mPmaxr\workspace.mat')
m = simparams.m;
n = simparams.n;

% Calc distance in initial positions x0_chaser and x0_target

trail_dist = vecnorm(x0_target(1:3) - x0_chaser(1:3)) * ndDist2km;
trail_time_hrs = time_past_chaser_target0 * ndTime2hrs;

x1 = reshape(x_opt, simparams.m, simparams.n);
traj1 = traj;
r_plot1 = traj.x_t(:,1:3);
simparams1 = simparams;





%% Traj 2 - 1 hour trail
% C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\20231005_0846.06_rdvz_nrho_2dv_1hrlead_det_60secMinIntTime
load('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\20231005_0846.06_rdvz_nrho_2dv_1hrlead_det_60secMinIntTime\workspace.mat')
trail_dist = vecnorm(x0_target(1:3) - x0_chaser(1:3)) * ndDist2km;
trail_time_hrs = time_past_chaser_target0 * ndTime2hrs;


x2 = reshape(x_opt, simparams.m, simparams.n);
traj2 = traj;
r_plot2 = traj.x_t(:,1:3);
simparams2 = simparams;


%% Traj 3 - 3 hour trail
% C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\20230928_1312.00_rdvz_nrho_2dv_s4EL_3hrLead_det_Qem8_6TcmMax_5seg_compareF
load('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\20230928_1312.00_rdvz_nrho_2dv_s4EL_3hrLead_det_Qem8_6TcmMax_5seg_compareF\workspace.mat')

trail_dist = vecnorm(x0_target(1:3) - x0_chaser(1:3)) * ndDist2km;
trail_time_hrs = time_past_chaser_target0 * ndTime2hrs;

x3 = reshape(x_opt, simparams.m, simparams.n);
traj3 = traj;
r_plot3 = traj.x_t(:,1:3);
simparams3 = simparams;


%% Traj 4 - 8 hour trail
% C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\20230928_1327.06_rdvz_nrho_2dv_s4EL_8hrLead_det_Qem8_6TcmMax_5seg_compareG
load('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\20230928_1327.06_rdvz_nrho_2dv_s4EL_8hrLead_det_Qem8_6TcmMax_5seg_compareG\workspace.mat')

trail_dist = vecnorm(x0_target(1:3) - x0_chaser(1:3)) * ndDist2km;
trail_time_hrs = time_past_chaser_target0 * ndTime2hrs;

x4 = reshape(x_opt, simparams.m, simparams.n);
traj4 = traj;
r_plot4 = traj.x_t(:,1:3);
simparams4 = simparams;


%% Traj 5 - 5 minute trail
%%%%% may need to create this one




%% Plot
figure;
x_init = simparams.x_init;
options = simparams.options;
x_target = simparams.x_target;


dynSys = simparams.dynSys;

if strcmp(dynSys,'2bp')
    [~, x_1] = ode113(@r2bp_de, [0,simparams.T0], x_init, options, mu);
    [~, x_2] = ode113(@r2bp_de, [0,simparams.T_target], x_target, options, mu);
elseif strcmp(dynSys,'cr3bp')
    [~, x_1] = ode113(@cr3bp_sFrame_nd_de, [0,simparams.T0], x_init, options, mu);
    [~, x_2] = ode113(@cr3bp_sFrame_nd_de, [0,simparams.T_target], x_target, options, mu);
end

plot3(x_1(:,1), x_1(:,2), x_1(:,3),'Color','Black')
plot3(x_2(:,1), x_2(:,2), x_2(:,3),'Color','Black')
hold on;
axis equal
grid on
plot3(r_plot1(:,1), r_plot1(:,2), r_plot1(:,3),'LineWidth',3);
plot3(r_plot2(:,1), r_plot2(:,2), r_plot2(:,3),'LineWidth',3);
plot3(r_plot3(:,1), r_plot3(:,2), r_plot3(:,3),'LineWidth',3);
plot3(r_plot4(:,1), r_plot4(:,2), r_plot4(:,3),'LineWidth',3);




plotNomDv(x1, traj1, simparams1)
plotNomDv(x2, traj2, simparams2)
plotNomDv(x3, traj3, simparams3)
plotNomDv(x4, traj4, simparams4)


