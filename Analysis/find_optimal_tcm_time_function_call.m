clc; close all; clear;

%% load a LEO coast to Hohmann xfer to GEO 
load('C:\Users\skell\OneDrive - USU\Documents\code_repos\multiseg_opt\sims\20220706_1627.55_leo2leo_smallP_objandconstraintGrad_testsnoSTT_EXAMPLE1\workspace.mat')
addpath(genpath('../'))
% Traj is in x_opt
x = x_opt;

%% Define P_initial
% Initial uncertainty
% Small
sig_pos = 10 / 1e3; % Position +/- 10 m in all 3 direction
sig_vel = .01 / 1e3; % Velocity +/- 1 cm/s in all 3 directions
% Medium
% sig_pos = 1000 / 1e3; % Position +/- 1 km in all 3 direction
% sig_vel = 1 / 1e3; % Velocity +/- 1 m/s in all 3 directions
% Large
% sig_pos = 10000 / 1e3; % Position +/- 10 km in all 3 direction
% sig_vel = 10 / 1e3; % Velocity +/- 10 m/s in all 3 directions
% Huge
% sig_pos = 100000 / 1e3; % Position +/- 100 km in all 3 direction
% sig_vel = 100 / 1e3; % Velocity +/- 100 m/s in all 3 directions


% sig_pos = 1; % Position +/- 10 m in all 3 direction
% sig_vel = 0; % Velocity +/- 1 cm/s in all 3 directions
% P_initial = diag([sig_pos^2 sig_pos^2 sig_pos^2 sig_vel^2 sig_vel^2 sig_vel^2]);
P_initial = diag([0 sig_pos^2 sig_pos^2  sig_vel^2 sig_vel^2 sig_vel^2]);

%% Call optimal_single_tcm function


[t_c_opt, min3sigmaDVC, corrGrad, dvC3sigma, t, t_zerox] = optimal_single_tcm(x, P_initial, mu);


figure
plot(t, corrGrad,'DisplayName','dTCM/dtc')
hold on
yline(0)
limit = trimmean(corrGrad,20);
ylim([-limit limit])
xlim([0 t(end)])
% 
plot(t_zerox, 0*t_zerox,'.','MarkerSize',25,'DisplayName','Zero crossings')
% 
plot(t_c_opt,0,'o','MarkerSize',10,'DisplayName','Zero corresponding to minimum')
% 
xlabel('Time (hrs)')
legend('Location','SouthEast')
title('Partial of 3\sigma TCM magnitude wrt TCM execution time')
grid on


figure
plot(t,dvC3sigma,'DisplayName','3\sigma TCM')
hold on
xlabel('Time (hrs)')
ylabel('Delta V (km/hr)')
limit = trimmean(dvC3sigma,20)*3;
ylim([0 limit])
xlim([0 t(end)])
plot(t_c_opt,min3sigmaDVC,'o','MarkerSize',10,'DisplayName','Minimum value (identified from plot 1)')
legend('Location','SouthEast')
grid on
title('3\sigma TCM magnitude')

figure
plot(t,dvC3sigma,'DisplayName','3\sigma TCM')
hold on
xlabel('Time (hrs)')
ylabel('Delta V (km/hr)')
limit = trimmean(dvC3sigma,20);
ylim([0 limit])
plot(t_c_opt,min3sigmaDVC,'o','MarkerSize',10,'DisplayName','Minimum value (identified from plot 1)')
legend('Location','SouthEast')
grid on
title('3\sigma TCM magnitude')




