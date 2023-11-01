clear; close all; clc;


%% Traj 1 - 12 hour trail
% C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\20231004_1516.52_rdvz_nrho_2dv_12hrlead_det_Qem8_4TcmMax_100mPmaxr

% load('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\20231004_1516.52_rdvz_nrho_2dv_12hrlead_det_Qem8_4TcmMax_100mPmaxr\workspace.mat')
load('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\20231031_1544.46_nrho_rdvz_det_startNearPerilune\workspace.mat')
colorblind = simparams.colorblind;

m = simparams.m;
n = simparams.n;

% Calc distance in initial positions x0_chaser and x0_target

trail_dist1 = vecnorm(x0_target(1:3) - x0_chaser(1:3)) * ndDist2km;
trail_time_hrs1 = time_past_chaser_target0 * ndTime2hrs;

plabel1 = strcat('Initial trail distance: ',32, num2str(round(trail_dist1)),' km; Initial delay: ',32,num2str(trail_time_hrs1),' hrs');

x1 = reshape(x_opt, simparams.m, simparams.n);
traj1 = traj;
r_plot1 = traj.x_t(:,1:3);
simparams1 = simparams;


sum(x1(7,2:end-1))*ndTime2days


%% Traj 2 - 1 hour trail
% C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\20231005_0846.06_rdvz_nrho_2dv_1hrlead_det_60secMinIntTime
load('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\20231005_0846.06_rdvz_nrho_2dv_1hrlead_det_60secMinIntTime\workspace.mat')
trail_dist2 = vecnorm(x0_target(1:3) - x0_chaser(1:3)) * ndDist2km;
trail_time_hrs2 = time_past_chaser_target0 * ndTime2hrs;

plabel2 = strcat('Initial trail distance: ',32, num2str(round(trail_dist2)),' km; Initial delay: ',32,num2str(trail_time_hrs2),' hrs');


x2 = reshape(x_opt, simparams.m, simparams.n);
traj2 = traj;
r_plot2 = traj.x_t(:,1:3);
simparams2 = simparams;

sum(x2(7,2:end-1))*ndTime2days
%% Traj 3 - 3 hour trail
% C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\20230928_1312.00_rdvz_nrho_2dv_s4EL_3hrLead_det_Qem8_6TcmMax_5seg_compareF
load('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\20230928_1312.00_rdvz_nrho_2dv_s4EL_3hrLead_det_Qem8_6TcmMax_5seg_compareF\workspace.mat')

trail_dist3 = vecnorm(x0_target(1:3) - x0_chaser(1:3)) * ndDist2km;
trail_time_hrs3 = time_past_chaser_target0 * ndTime2hrs;

plabel3 = strcat('Initial trail distance: ',32, num2str(round(trail_dist3)),' km; Initial delay: ',32,num2str(trail_time_hrs3),' hrs');


x3 = reshape(x_opt, simparams.m, simparams.n);
traj3 = traj;
r_plot3 = traj.x_t(:,1:3);
simparams3 = simparams;
sum(x3(7,2:end-1))*ndTime2days

%% Traj 4 - 8 hour trail
% C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\20230928_1327.06_rdvz_nrho_2dv_s4EL_8hrLead_det_Qem8_6TcmMax_5seg_compareG
load('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\20230928_1327.06_rdvz_nrho_2dv_s4EL_8hrLead_det_Qem8_6TcmMax_5seg_compareG\workspace.mat')

trail_dist4 = vecnorm(x0_target(1:3) - x0_chaser(1:3)) * ndDist2km;
trail_time_hrs4 = time_past_chaser_target0 * ndTime2hrs;

plabel4 = strcat('Initial trail distance: ',32, num2str(round(trail_dist4)),' km; Initial delay: ',32,num2str(trail_time_hrs4),' hrs');


x4 = reshape(x_opt, simparams.m, simparams.n);
traj4 = traj;
r_plot4 = traj.x_t(:,1:3);
simparams4 = simparams;
sum(x4(7,2:end-1))*ndTime2days

%% Traj 5 - 5 minute trail
% C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\20231031_1508.43_nrho_rdvz_det_startNearPerilune
load('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\20231031_1508.43_nrho_rdvz_det_startNearPerilune\workspace.mat')

trail_dist5 = vecnorm(x0_target(1:3) - x0_chaser(1:3)) * ndDist2km;
trail_time_hrs5 = time_past_chaser_target0 * ndTime2hrs*60;

plabel5 = strcat('Initial trail distance: ',32, num2str(round(trail_dist5)),' km; Initial delay: ',32,num2str(trail_time_hrs5),' mins');


x5 = reshape(x_opt, simparams.m, simparams.n);
traj5 = traj;
r_plot5 = traj.x_t(:,1:3);
simparams5 = simparams;
sum(x5(7,2:end-1))*ndTime2days
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

% plot3(x_1(:,1), x_1(:,2), x_1(:,3),'Color','Black')
plot3(x_2(:,1), x_2(:,2), x_2(:,3),'Color','Black','DisplayName','9:2 NRHO')
hold on;







axis equal
grid on

plot3(r_plot1(:,1), r_plot1(:,2), r_plot1(:,3),'LineWidth',3,'Color',colorblind(1,:),'DisplayName',plabel1);
plot3(r_plot4(:,1), r_plot4(:,2), r_plot4(:,3),'LineWidth',3,'Color',colorblind(4,:),'DisplayName',plabel4);

plot3(r_plot3(:,1), r_plot3(:,2), r_plot3(:,3),'LineWidth',3,'Color',colorblind(3,:),'DisplayName',plabel3);
plot3(r_plot2(:,1), r_plot2(:,2), r_plot2(:,3),'LineWidth',3,'Color',colorblind(2,:),'DisplayName',plabel2);

plot3(r_plot5(:,1), r_plot5(:,2), r_plot5(:,3),'LineWidth',3,'Color',colorblind(5,:),'DisplayName',plabel5);
legend();




plotNomDv(x1, traj1, simparams1)
plotNomDv(x2, traj2, simparams2)
plotNomDv(x3, traj3, simparams3)
plotNomDv(x4, traj4, simparams4)
plotNomDv(x5, traj5, simparams5)
    xlabel('X (nd)')
    ylabel('Y (nd)')
    zlabel('Z (nd)')
        view([80, 12]); %leo-nri view

%         view([90 0])
%         view([0 0])


%% Subplot figure
close all;


solfig = figure;


s1=subplot(1,2,1)

% plot3(x_1(:,1), x_1(:,2), x_1(:,3),'Color','Black')
plot3(x_2(:,1), x_2(:,2), x_2(:,3),'Color','Black','DisplayName','9:2 NRHO')
hold on;







axis equal
grid on

plot3(r_plot1(:,1), r_plot1(:,2), r_plot1(:,3),'LineWidth',3,'Color',colorblind(1,:),'DisplayName',plabel1);
plot3(r_plot4(:,1), r_plot4(:,2), r_plot4(:,3),'LineWidth',3,'Color',colorblind(4,:),'DisplayName',plabel4);

plot3(r_plot3(:,1), r_plot3(:,2), r_plot3(:,3),'LineWidth',3,'Color',colorblind(3,:),'DisplayName',plabel3);
plot3(r_plot2(:,1), r_plot2(:,2), r_plot2(:,3),'LineWidth',3,'Color',colorblind(2,:),'DisplayName',plabel2);

plot3(r_plot5(:,1), r_plot5(:,2), r_plot5(:,3),'LineWidth',3,'Color',colorblind(5,:),'DisplayName',plabel5);
% legend();




plotNomDv(x1, traj1, simparams1)
plotNomDv(x2, traj2, simparams2)
plotNomDv(x3, traj3, simparams3)
plotNomDv(x4, traj4, simparams4)
plotNomDv(x5, traj5, simparams5)
    xlabel('X (nd)')
    ylabel('Y (nd)')
    zlabel('Z (nd)')
%         view([80, 12]); %leo-nri view

        view([90 0])
%         view([0 0])

hold off
s2=subplot(1,2,2)

% plot3(x_1(:,1), x_1(:,2), x_1(:,3),'Color','Black')
plot3(x_2(:,1), x_2(:,2), x_2(:,3),'Color','Black','DisplayName','9:2 NRHO')
hold on;







axis equal
grid on

plot3(r_plot1(:,1), r_plot1(:,2), r_plot1(:,3),'LineWidth',3,'Color',colorblind(1,:),'DisplayName',plabel1);
plot3(r_plot4(:,1), r_plot4(:,2), r_plot4(:,3),'LineWidth',3,'Color',colorblind(4,:),'DisplayName',plabel4);

plot3(r_plot3(:,1), r_plot3(:,2), r_plot3(:,3),'LineWidth',3,'Color',colorblind(3,:),'DisplayName',plabel3);
plot3(r_plot2(:,1), r_plot2(:,2), r_plot2(:,3),'LineWidth',3,'Color',colorblind(2,:),'DisplayName',plabel2);

plot3(r_plot5(:,1), r_plot5(:,2), r_plot5(:,3),'LineWidth',3,'Color',colorblind(5,:),'DisplayName',plabel5);
lg = legend('Location','southoutside');




plotNomDv(x1, traj1, simparams1)
plotNomDv(x2, traj2, simparams2)
plotNomDv(x3, traj3, simparams3)
plotNomDv(x4, traj4, simparams4)
plotNomDv(x5, traj5, simparams5)
    xlabel('X (nd)')
    ylabel('Y (nd)')
    zlabel('Z (nd)')
%         view([80, 12]); %leo-nri view

%         view([90 0])
        view([0 0])


fig = gcf;
lg.Position(2)=.025
lg.Position(1)=.23



s1.OuterPosition(4)=.75
s1.Position(2) = .35

s2.OuterPosition(4)=.75
s2.Position(2) = .35

solfig.CurrentAxes.Title=text(1,1,'nrho_det_opt_trajs');
solfig.CurrentAxes.Title.Visible="off";


outputPath = 'C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\Analysis\nrho_rdvz';
saveallfigs(outputPath,0)

