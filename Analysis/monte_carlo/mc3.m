% try and incorporate dispersions only in the LVLH directions to understand
% the impact
cd('C:\Users\skell\OneDrive - USU\Documents\code_repos\multiseg_opt_manual_tcm\Analysis\monte_carlo')

clear;
clear global;
close all; 
clc;
format longg; 
addpath(genpath('../../'));


% load('C:\Users\skell\OneDrive - USU\Documents\code_repos\multiseg_opt_manual_tcm\sims\20221103_1138.29_leo45_to_leoSunSync_testing_deterministic_incChangeAllDV1_targetXFN_freeTCMloc\workspace.mat')
% load('C:\Users\skell\OneDrive - USU\Documents\code_repos\multiseg_opt_manual_tcm\sims\20221101_1623.28_leo45_to_leoSunSync_robust_10kmR_incChangeAllDV1_J\workspace.mat')
% load('C:\Users\skell\OneDrive - USU\Documents\code_repos\multiseg_opt_manual_tcm\sims\20221101_1333.01_leo28_to_geo0_J\workspace.mat')
% load('C:\Users\skell\OneDrive - USU\Documents\code_repos\multiseg_opt_manual_tcm\sims\20221108_1442.58_leo45toleo45_ROBUST\workspace.mat')
% load('C:\Users\skell\OneDrive - USU\Documents\code_repos\multiseg_opt_manual_tcm\sims\20221121_2001.27_leo45toleo45_dv1tcm1_concurrent\workspace.mat')
% load('C:\Users\skell\OneDrive - USU\Documents\code_repos\multiseg_opt_manual_tcm\sims\20221122_1346.21_leo45toleo45_dv1tcm1_concurrent_newStatsiDv1\workspace.mat')
% load('C:\Users\skell\OneDrive - USU\Documents\code_repos\multiseg_opt_manual_tcm\sims\20221121_1723.15_leo45toleo45_final_dvCTcm_unitVec\workspace.mat')
% load('C:\Users\skell\OneDrive - USU\Documents\code_repos\multiseg_opt_manual_tcm\sims\20221122_1755.47_leo28_to_geo0_beforeZeroStart_robust\workspace.mat')
% load('C:\Users\skell\OneDrive - USU\Documents\code_repos\multiseg_opt_manual_tcm\sims\20221122_1805.18_leo28_to_geo0_beforeZeroStart_dv1tcm1concurrent\workspace.mat')
% load('C:\Users\skell\OneDrive - USU\Documents\code_repos\multiseg_opt_manual_tcm\sims\20221122_1815.36_leo28_to_geo0_beforeZeroStart_dv1tcm1concurrent_iDV1stats\workspace.mat')
% load('C:\Users\skell\OneDrive - USU\Documents\code_repos\multiseg_opt_manual_tcm\sims\20221122_1844.51_leo28_to_geo0_preZeroStart_dv1tcm1concurrent_trTcmStats\workspace.mat')
% load('C:\Users\skell\OneDrive - USU\Documents\code_repos\multiseg_opt_manual_tcm\sims\20221122_1902.09_leoCrit_to_leoSunSync_dv1tcm1Tied_allplanechangeDV1\workspace.mat')
% load('C:\Users\skell\OneDrive - USU\Documents\code_repos\multiseg_opt_manual_tcm\sims\20221122_1905.50_leoCrit_to_leoSunSync_dv1tcm1Tied_allplanechangeDV1_idvStats\workspace.mat')
%  load('C:\Users\skell\OneDrive - USU\Documents\code_repos\multiseg_opt_manual_tcm\sims\20221123_0927.48_leoCrit_to_leoSunSync_dv1tcm1Tied_allplanechangeDV1_idvR_trV\workspace.mat')
 load('C:\Users\skell\OneDrive - USU\Documents\code_repos\multiseg_opt_manual_tcm\sims\20221123_1345.48_leo28_to_geo0_robust\workspace.mat')
 
% simparams.sig_pos = 100;
% simparams.sig_vel = .1 * 3600;
% simparams.P_initial = diag([simparams.sig_pos^2 simparams.sig_pos^2 simparams.sig_pos^2 simparams.sig_vel^2 simparams.sig_vel^2 simparams.sig_vel^2]);



sigr_in_track = 10;
sigr_out_of_plane = 10;
sigr_altitude = 10;
sig_r = [sigr_in_track; sigr_out_of_plane; sigr_altitude];


sigv_in_track = .01 * 3600;
sigv_out_of_plane = .01 * 3600;
sigv_altitude = .01*3600;
sig_v = [sigv_in_track; sigv_out_of_plane; sigv_altitude];


r0 = x_opt(1:3,1);
v0 = x_opt(4:6,1);

% out of plane
j_hat = cross(r0,v0) / norm( cross(r0,v0));

% altitude
k_hat = r0 / norm(r0);

% in track
% i_hat = v0 / norm(v0); % only for circular
i_hat = cross(j_hat, k_hat);

% T_i_lvlh = [i_hat'; j_hat'; k_hat'];
% T_lvlh_i = T_i_lvlh';

T_lvlh_i = [i_hat'; j_hat'; k_hat'];
T_i_lvlh = T_lvlh_i';

bdT = zeros(6,6);
bdT(1:3,1:3) = T_i_lvlh;
bdT(4:6,4:6) = T_i_lvlh;
simparams.P_initial = bdT * diag([sig_r.^2; sig_v.^2]) * bdT';



%% Monte carlo simulation of trajectory x_opt

x_opt = reshape(x_opt,simparams.m, simparams.n);
% Simulate trajectory and determine miss distance at the target




% Propagate unperturbed initial trajectory, calculate dv1 and dv2

[stm_i, x_i_f, x_t, stm_t, t] = createStateStmHistory(x_opt, simparams);
[~, deltaVs_nom] = calcDeltaV(x_opt,x_i_f,simparams);

maneuver_times = [];
for i = 1:length(simparams.maneuverSegments)
    maneuver_times(i) = sum(x_opt(7,1:simparams.maneuverSegments(i)-1));
end

% Determine target states and times (maneuver/final trajectory...I suppose
% tcm_target is the only one that matters)
final_target = simparams.x_target;
final_target_time = sum(x_opt(7,:));

% simparams.target_final_maneuver = 1; %%%% CHANGED TO TARGET FINAL MANEUVER

if simparams.target_final_maneuver
    tcm_target = x_opt(1:6,simparams.maneuverSegments(end)) - [zeros(3,1); deltaVs_nom(:,2)];
    tcm_target_time = sum(x_opt(7,1:simparams.maneuverSegments(end)-1));
else
    tcm_target = final_target;
    tcm_target_time = final_target_time;
end

% Get the tcm execution time
[tcm_min, tcm_time, dvR3sigma_tr, dvV3sigma_tr, dvR3sigma_i, dvV3sigma_i, P_tcm1, P_tcm2] = tcmPair_rv(x_opt, t, stm_t, deltaVs_nom, simparams);

%%%%%% CHANGED TO ALIGN TCM WITH FIRST DELTA V
% tcm_time = maneuver_times(1);


fig_inertial = figure;
plot3(x_t(:,1),x_t(:,2),x_t(:,3),'LineWidth',2,'Color','Black')
axis equal
hold on
grid on
plot3(tcm_target(1),tcm_target(2),tcm_target(3),'.','markersize',25)


% fig_components = figure;
% hold on
% axis equal
% grid on


tcm_norm_sav = [];
tcm_r_el_sav = [];
tcm_v_el_sav = [];


%% Start random variable part
for j = 1:100

    
    % Synthesize random initial state realization
%     r0_err = randn(3,1).*[simparams.sig_pos];
    r0_err = randn .* (T_i_lvlh * sig_r);

%     v0_err = randn(3,1)*simparams.sig_vel;
    v0_err = randn .* (T_i_lvlh * sig_v);
    
%     r0_err = zeros(3,1);
%     v0_err = zeros(3,1);
    del_x_0 = [r0_err; v0_err];
    x0 = x_opt(1:6,1) + del_x_0;

    tcm1_dv = zeros(3,1);

    correction_needed = 1;
    cnt = 0;

    while correction_needed
        
        % Propagate initial error state thru trajectory with no TCM (tcm
        % magnitudes zero)
        event_times = [maneuver_times,tcm_time, tcm_target_time];
        event_dvs = [deltaVs_nom, tcm1_dv, zeros(3,1)];
        [x_tt, stm_tt, tt, x_currt] = mcProp(x0,event_times,event_dvs,simparams);
        
        
        
        % Use STM to correct del_r via tcm_dv executed at tcm_time
        tcm_idx(1) = find( tcm_time == tt );
        tcm_idx(2) = find( tcm_target_time == tt );
    
        % Determine tcm target dispersion
        del_x_n_num = x_tt(:,end) - tcm_target;
        
        stmN0 = stm_tt(:,:,tcm_idx(2));
        J = [zeros(3,3), eye(3,3); -eye(3,3), zeros(3,3)];
        stmC0 = stm_tt(:,:,tcm_idx(1));
        stm0C = -J * stmC0' * J; % symplectic inverse
        stmNC = stmN0 * stm0C;
        
        T = [-inv( stmNC(1:3,4:6) ) * stmNC(1:3,1:3), -eye(3)];
        
    
        del_x_c = x_tt(:,tcm_idx(1)) - x_t(t==tcm_time,:)';
        tcm1_dv = T * del_x_c;
    
        if simparams.nom_dvctied
            tcm1_norm_sav(j) = norm(event_dvs(:,1) - tcm1_dv) - norm(event_dvs(:,1));
            tcm_norm_sav(j) = tcm1_norm_sav(j);
        else
            tcm1_norm_sav(j) = norm(tcm1_dv);
            tcm_norm_sav(j) = tcm1_norm_sav(j);
            tcm_r_el_sav(1:3,j) = tcm1_dv;
        end

        
        % test out stm linear propagation of del_x_0 to del_x_f
        del_x_n_lin = stmNC * stmC0 * del_x_0;
    %     linPropErrorPercent = (del_x_n_num - del_x_n_lin)./del_x_n_num * 100
        
        
        
        %% Propagate again with correction, see how close it comes
        event_times = [maneuver_times, tcm_time, tcm_target_time];
        event_dvs = [deltaVs_nom, tcm1_dv, zeros(3,1)];
        [x_te, stm_te, te, x_curre] = mcProp(x0,event_times,event_dvs,simparams);
        
        % tcm target was calculated without the nominal deltaV 2
%         del_xe = tcm_target - x_curre;
%         del_xe = x_curre - x_opt(1:6,simparams.maneuverSegments(end));

        del_xe = x_opt(1:6,simparams.maneuverSegments(end)) - x_curre;
        tcm2_dv = del_xe(4:6);
        
%         del_xe_f = stmNC * (del_x_0 + [zeros(3,1); tcm_dv]) %%%% error if tcm does not happen at the start of the traj
    
%         tcm_norm_sav(j) = tcm_norm_sav(j) + norm(del_xe(4:6));
        % Changing this part to the difference from the sum of the tcm and
        % the nominal minus the nominal, vs taking the tcm norm separately
%         tcm_norm_sav(j) = tcm_norm_sav(j) + norm(del_xe(4:6));


%         tcm2_dv = deltaVs_nom(:,2) - del_xe(4:6);
%         diff = norm(tcm2_dv) - norm(deltaVs_nom(:,2));
%         tcm2_norm_sav(j) = diff;
%         tcm_norm_sav(j) = tcm_norm_sav(j) + diff;



        tcm2_norm_sav(j) = norm(deltaVs_nom(:,2) + tcm2_dv) - norm(deltaVs_nom(:,2));
        tcm_norm_sav(j) = tcm_norm_sav(j) + tcm2_norm_sav(j);
        
        
        % the velocity correction TCM is really a correction to DV2:
        


        if norm(del_xe(1:3)) < .1
            norm(del_xe(1:3))
            correction_needed = 0;
            
        end
        cnt = cnt+1;

        if cnt > 3
            tcm1_dv_new = fsolve( @(tcm1_dv)fsolve_tcm_target_function(tcm1_dv, simparams, x_opt, deltaVs_nom, maneuver_times, tcm_time, tcm_target_time, x0), tcm1_dv );
            tcm1_dv = tcm1_dv_new;

            if simparams.nom_dvctied
                tcm1_norm_sav(j) = norm(event_dvs(:,1) - tcm1_dv) - norm(event_dvs(:,1));
                tcm_norm_sav(j) = tcm1_norm_sav(j);
            else
                tcm1_norm_sav(j) = norm(tcm1_dv);
                tcm_norm_sav(j) = tcm1_norm_sav(j);
                tcm_r_el_sav(1:3,j) = tcm1_dv;
            end


            event_times = [maneuver_times, tcm_time, tcm_target_time];
            event_dvs = [deltaVs_nom, tcm1_dv_new, zeros(3,1)];
            [x_te, stm_te, te, x_curre] = mcProp(x0,event_times,event_dvs,simparams);

            del_xe = x_opt(1:6,simparams.maneuverSegments(end)) - x_curre;
            tcm2_dv = del_xe(4:6);

            tcm2_norm_sav(j) = norm(deltaVs_nom(:,2) + tcm2_dv) - norm(deltaVs_nom(:,2));
            tcm_norm_sav(j) = tcm_norm_sav(j) + tcm2_norm_sav(j);

            correction_needed = 0;
        end


    % Propagate one final time to the end of the trajectory with both
    % corrections included
   
    if simparams.nom_dvctied
        event_times = [maneuver_times, sum(x_opt(7,:))];
        event_dvs = [deltaVs_nom(:,1) + tcm1_dv, deltaVs_nom(:,2) + tcm2_dv, zeros(3,1)];
    else
        event_times = [maneuver_times, tcm_time, tcm_target_time, sum(x_opt(7,:))];
        event_dvs = [deltaVs_nom(:,1), deltaVs_nom(:,2) + tcm2_dv, tcm1_dv, zeros(3,2)];
    end

    [x_tf, stm_tf, tf, x_currf] = mcProp(x0,event_times,event_dvs,simparams);

%     simparams.x_target - x_currf
%     diff
        
    end




    %% plot
    
%     plot3(x_tt(1,:),x_tt(2,:),x_tt(3,:))
    
    figure(fig_inertial)
    plot3(x_tf(1,:),x_tf(2,:),x_tf(3,:))
    
%     figure(fig_components)
%     plot(te, x_te(1,:),'Color','Red')
%     plot(te, x_te(2,:),'Color','Black')
%     plot(te, x_te(3,:),'Color','Cyan')
   

end

mean_tcm = mean(tcm_norm_sav)

threesigma_tcm = 3 * std(tcm_norm_sav)
tcm_min
tcm_min / (threesigma_tcm + mean_tcm)






figure
histogram(tcm1_norm_sav,'Normalization','probability')
title('TCM1 histogram')
hold on
xline(dvR3sigma_tr,'LineWidth',2,'Color','Red')
xline(dvR3sigma_i,'LineWidth',2,'Color','Black')
legend('','Trace 3\sigma','i_d_v method')

figure
histogram(tcm2_norm_sav,'Normalization','probability')
title('TCM2 histogram')
hold on
xline(dvV3sigma_tr,'LineWidth',2,'Color','Red')
xline(dvV3sigma_i,'LineWidth',2,'Color','Black')
legend('','Trace 3\sigma','i_\Delta_V method')


if ~simparams.nom_dvctied
    testCovR = cov(tcm_r_el_sav');
    tcm_r_tr = sqrt(trace(testCovR));
end
% 
% testCovV = cov(tcm_v_el_sav');
% tcm_v_tr = sqrt(trace(testCovV));

% testCovBoth = cov(tcm_r_el_sav'+tcm_v_el_sav');
% tcm3sigtotal = 3*sqrt(trace(testCovBoth));



%%%%% The difference is pretty large



cd('C:\Users\skell\OneDrive - USU\Documents\code_repos\multiseg_opt_manual_tcm')
saveallfigs(outputPath,0)



sigma_i = sqrt(eig(P_tcm2));

%% chi squared testing
% xi = 0:100:100000;
% f = zeros(1,length(xi));
% for i = 1:length(sigma_i)
%     prod = 1;
%     for j = 1:length(sigma_i)
%         if j ~= i
%             prod = prod * (1 - sigma_i(j)^2/sigma_i(i)^2);
%         end
%     end
%     f = f + exp(-xi/sigma_i(i)^2)./(sigma_i(i)^2 * prod);
% end

% figure
% plot(xi,f)




