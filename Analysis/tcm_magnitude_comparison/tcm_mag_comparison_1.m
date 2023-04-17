

% cd('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\Analysis\tcm_magnitude_comparison');
clear;
clear global;
close all; 
clc;
format longg; 



% load('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\sims\from_aries\20230404_2216.25_cr3bp_leo_mlo_3dv_midcourseTcmTarget_robust\workspace.mat');
load('workspace.mat');


% addpath(genpath('../../'));
% addpath(genpath('../robust/'));



simparams.R_dv = simparams.R;
[tcm_time,tcm_idx,min_tcm_dv, P_i_minus, P_i_plus] = opt_multiple_tcm(x_opt, t, t_s, stm_t, simparams);


P0 = P_i_plus(:,:,3);
% P0 = simparams.P_initial;


x = reshape(x_opt,simparams.m,simparams.n);
x = x(:,simparams.P_constrained_nodes(1):simparams.P_constrained_nodes(2)-1);
simparams.n = size(x,2);


[stm_i, stt_i, x_i_f, x_t, stm_t, stt_t_i, t, t_s, stm_t_i]  = createStateStmSttHistory(x, simparams);


stmN0 = stm_t(:,:,end);

% tcm_rss = zeros(1,length(t)-1);
% tcm_p44 = zeros(1,length(t)-1);
% tcm_mc = zeros(1,length(t)-1);

n = 1000;
x0  = x(1:6,1);
x_target = x_t(end,:)';

range = 1:25:length(t)-1;
length_range = length(range);

tcm_rss = zeros(1,length_range);
tcm_p44 = zeros(1,length_range);
tcm_mc = zeros(1,length_range);

% range = 1:100:length(t)-1;
% for i = 1:100:length(t)-1
parfor j = 1:length_range
% for j = 12
    j
    i = range(j);

    stmC0 = stm_t(:,:,i);
    stm0C = invert_stm(stmC0, simparams);

    stmNC = stmN0 * stm0C;

    P = stmC0 * P0 * stmC0';

    T = [-inv( stmNC(1:3,4:6) ) * stmNC(1:3,1:3), -eye(3)];

    P_tcm = T * P * T' + simparams.R;

    tcm_rss(j) = 3*sqrt(trace(P_tcm));
    tcm_p44(j) = DVmag99p7(P_tcm);

    tcm_time = t(i);
    tcm_mc(j) = monteCarloTcmMagnitude(x0, x_target, t, tcm_time, stm_t, P0, n, 99.7, simparams);

end





[min_rss, min_rss_idx] = min(tcm_rss);
[min_p44, min_p44_idx] = min(tcm_p44);
[min_mc, min_mc_idx] = min(tcm_mc);






% tcm_time = t(min_p44_idx);


% tcm_dv_mag = monteCarloTcmMagnitude(x0, x_target, t, tcm_time, stm_t, P0, n, 99.7, simparams)


plot_logical = tcm_rss > 0;



figure
% plot(t(1:end-1),tcm_rss,'LineWidth',2,'DisplayName','RSS');
% plot(t(plot_logical),tcm_rss(plot_logical),'LineWidth',2,'DisplayName','RSS');
plot(t(range),tcm_rss,'LineWidth',2,'DisplayName','RSS');
hold on;
% plot(t(1:end-1),tcm_p44,'LineWidth',2,'DisplayName','P44');
% plot(t(plot_logical),tcm_p44(plot_logical),'LineWidth',2,'DisplayName','P44');
plot(t(range),tcm_p44,'LineWidth',2,'DisplayName','P44');

% plot(t(1:end-1),tcm_mc,'LineWidth',2,'DisplayName','MC');
% plot(t(plot_logical),tcm_mc(plot_logical),'LineWidth',2,'DisplayName','MC');
plot(t(range),tcm_mc,'LineWidth',2,'DisplayName','MC');

% plot(t(min_rss_idx), min_rss,'.','MarkerSize',25,'DisplayName','Min RSS')
% 
% plot(t(min_p44_idx), min_p44,'.','MarkerSize',25,'DisplayName','Min P44')
% 
% plot(t(min_mc_idx), min_mc,'.','MarkerSize',25,'DisplayName','Min MC')


ylim([0 1]);
ylabel('TCM magnitude (ND velocity)')
xlabel('TCM execution time along trajectory (ND time)')
legend()


% figure
% plotMultiSegTraj(x, x_t, t_s, simparams);





%% Save

formatOut = 'yyyymmdd_HHMM.SS';
dateString = datestr(now,formatOut);
outputPath = strcat('./sims/',dateString);
mkdir(outputPath);
% Save workspace    
save(strcat('./',outputPath,'/workspace.mat'));

saveallfigs(strcat('./',outputPath),0)






















