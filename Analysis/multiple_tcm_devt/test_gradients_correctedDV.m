%% Multiple TCM stuff
% reference derivation in lyx document: stm_target_vel_correction

clear;
clear global;
% close all; 
clc;
format longg;
cd('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust');
addpath(genpath('./'));

%%
% init_fn = 'C:\Users\skell\OneDrive - USU\Documents\code_repos\robust_2body_tcmerror\init_traj_files\initialize_simulation_parameters_leo28_to_geo0';
% init_fn = './init_traj_files/init_simparams_cr3bp_leo_to_mlo_3dv';
init_fn = './init_traj_files/init_simparams_cr3bp_leoinclined_lloflyby_dro3_x0detOpt_3dv';



run(init_fn);
% simparams = generateInitialXfer(simparams);
%% ONLY If using a previous reference trajectory as the initial guess:
% load('.\init_traj_files\initial_guesses\tli_llo_deterministic_opt.mat')
% load('.\init_traj_files\initial_guesses\polar_llo_to_nrho_apolune.mat')
% simparams.x0 = x_opt;

%%

x = simparams.x0;
dx = sqrt(eps);

[stm_i, stt_i, x_i_f, x_t, stm_t, stt_t_i, t, t_s, stm_t_i]  = createStateStmSttHistory(x, simparams);
[deltaV, deltaVs_nom] = calcDeltaV(simparams.x0,x_i_f,simparams);


% % Using saved dynamics, calculate total impulsive delta V
% [deltaV, deltaVs_nom] = calcDeltaV(x, x_i_f, simparams);

% Using saved dynamics, calculate the TCMs
[tcm_time, tcm_idx, min_tcm_dv, P_i_minus, P_i_plus, tcm_dv_each] = opt_multiple_tcm(x, deltaVs_nom, t, t_s, stm_t, stm_t_i, simparams);

% 
% [~, ~, tcm_dv_each, P_i_minus, P_i_plus] = calc_covariance_tcmdv(x, t, stm_t, tcm_time, simparams);

%%
R = simparams.R;
G = [zeros(3,3); eye(3,3)];

x = reshape(x,simparams.m,simparams.n);




%% Testing the gradient of a specific TCM
k=1; %%%%%% which TCM we're testing from tcm_time / tcm_idx

range(1) = 1;
second_maneuver_seg = simparams.maneuverSegments(2);
% range(2) = find(t_s==second_maneuver_seg,1);
range(2) = length(t);

[event_times, event_indicator] = define_events_v2(x(:), t, tcm_time, simparams);
event_idx_logical = logical(sum(t'==event_times', 1));    
event_idxs = find(event_idx_logical);

% %%%%%% FOR TESTING %%%%%%
% % Forcing a tcm time to be the same as a nominal maneuver 
% tcm_time(2) = event_times(3);
% % [event_times, event_is_tcm] = define_events(x(:), t, tcm_time, simparams);
% 
% [event_times, event_indicator] = define_events_v2(x(:), t, tcm_time, simparams);
% % 
% % 
% [tcm_timeo, ~, min_tcm_dv, P_i_minus, P_i_plus, tcm_dv_each] = opt_multiple_tcm_fdGradientTesting(x, t, t_s, stm_t, tcm_time, simparams);
% 
% %%%%%% END TESTING %%%%%%


num_events = length(event_indicator);
% event_times = zeros(num_events,1);
% event_times(event_is_tcm) = tcm_time;
% event_time_not_tcm = [];
% 
% for i = 1:sum(~event_is_tcm)
%     event_time_not_tcm(end+1) = sum(x(7,1:simparams.P_constrained_nodes(i) - 1));
% end
% 
% event_times(~event_is_tcm) = event_time_not_tcm;

 
% Logically indexing - return the event times that are not TCMs that are
% after the TCM
target_times = event_times( tcm_time(k) < event_times & event_indicator ~=1 );
target_time = target_times(1);
target_idx = find(t == target_time);




% stmNCk = dynCellCombine(t, t_s, tcm_idx(k), target_idx, simparams, stm_t_i);
stmNCk = dynCellCombine(t, t_s, event_idxs(k), target_idx, simparams, stm_t_i);


Tk = [-inv( stmNCk(1:3,4:6) ) * stmNCk(1:3,1:3), -eye(3)];


Pck_minus = P_i_minus(:,:,event_times==tcm_time(k));
sig_dvk = sqrt(trace(Tk * Pck_minus * Tk'));

Ptcm = Tk * Pck_minus * Tk';


tcmR_dv_total = sum(tcm_dv_each(1:end-1));
tcm1_dv = tcm_dv_each(1);
tcm2_dv = tcm_dv_each(2);
tcm3_dv = tcm_dv_each(3);
% tcm4_dv = tcm_dv_each(4);
% tcm5_dv = tcm_dv_each(5);

dtcmR_fd = zeros(length(x(:)),1);
dtcmV_fd = zeros(length(x(:)),1);
dtcm_total_fd=zeros(length(x(:)),1);
dTkdxi_fd = zeros(3,6,length(x(:)));
dstmNCk_fd = zeros(6,6,length(x(:)));


nomDVindex = 1;
DV = deltaVs_nom(:,nomDVindex);
DV_mag = vecnorm(DV);
i_DV = DV/DV_mag;

sig_corrected_DV = i_DV' * Ptcm * i_DV;




% t20_idx = find(t==x(7,1));
% t2f_idx = find(t==sum(x(7,1:2)));
% [stm_20_tc2, stt_20_tc2]  = dynCellCombine(t, t_s, tcm_idx(k), t20_idx, simparams, stm_t_i, stt_t_i);
% stm_2f_20  = dynCellCombine(t, t_s, t20_idx, t2f_idx, simparams, stm_t_i);



% if k == 1
%     tClast_idx = 1;
% else
%     tClast_idx = find(t==tcm_time(k-1));
% end
% tCk_idx = find(t==tcm_time(k));
% stmCkClast = dynCellCombine(t, t_s, tClast_idx, tCk_idx, simparams, stm_t_i);



parfor i = 1:simparams.m*simparams.n
% for i = 1:simparams.m*simparams.n
    i
    x1dx = x(:);
    x1dx(i) = x(i) + dx;
    [stm_i1dx, stt_i1dx, x_i_f1dx, x_t1dx, stm_t1dx, stt_t_i1dx, t1dx, t_s1dx, stm_t_i1dx]  = createStateStmSttHistory(x1dx, simparams);
%     [x_t1dx, stm_t1dx, stt_t_i1dx, t1dx, t_s1dx] = addToStateStmSttHistory(x_t1dx, stm_t1dx, stt_t_i1dx, t1dx, t_s1dx, [tcm_time, target_time], simparams);
    [x_t1dx, stm_t1dx, stm_t_i1dx stt_t_i1dx, t1dx, t_s1dx] = addToStateStmSttHistory(x_t1dx, stm_t1dx, stm_t_i1dx, stt_t_i1dx, t1dx, t_s1dx, [event_times], simparams);

    [deltaV1dx, deltaVs_nom1dx] = calcDeltaV(x1dx,x_i_f1dx,simparams);
%     [tcm_time1dx, ~, min_tcm_dv1dx, P_i_minus1dx, P_i_plus1dx, tcm_dv_each1dx] = opt_multiple_tcm_v2_fdGradientTesting(x1dx, deltaVs_nom1dx, t1dx, t_s1dx, stm_t1dx, stm_t_i1dx, tcm_time, simparams);

    [~, min_tcm_dv1dx, tcm_dv_each1dx, P_i_minus1dx, P_i_plus1dx] = calc_covariance_tcmdv_v3(x1dx(:), t1dx, t_s1dx, stm_t1dx, stm_t_i1dx, tcm_time, 1, deltaVs_nom1dx, simparams.P_initial, simparams);


%     [tcm_time1dx,tcm_idx1dx,min_tcm_dv1dx,~,~,tcm_num_option_DVs1dx] = opt_multiple_tcm(x1dx, t1dx, stm_t1dx, simparams);
%     [tcm_time,tcm_idx,min_tcm_dv,tcm_time_cell,tcm_idx_cell,tcm_num_option_DVs] 
%     [~, tcm_dv_total1dx, tcm_dv1dx] = calc_covariance_tcmdv(x1dx, t1dx, stm_t1dx, tcm_time, simparams);
%     [~, ~, tcm_dv_each1dx, P_i_minus1dx, P_i_plus1dx] = calc_covariance_tcmdv(x1dx, t1dx, stm_t1dx, tcm_time, simparams);

    tcmR_dv_total1dx = sum(tcm_dv_each1dx(1:end-1));

    dtcm1_fd(i) = (tcm_dv_each1dx(1) - tcm_dv_each(1))/dx;
    dtcm2_fd(i) = (tcm_dv_each1dx(2) - tcm_dv_each(2))/dx;
    dtcm3_fd(i) = (tcm_dv_each1dx(3) - tcm_dv_each(3))/dx;
    dtcm4_fd(i) = (tcm_dv_each1dx(4) - tcm_dv_each(4))/dx;
    dtcm5_fd(i) = (tcm_dv_each1dx(5) - tcm_dv_each(5))/dx;
    dtcm6_fd(i) = (tcm_dv_each1dx(6) - tcm_dv_each(6))/dx;
    dtcm7_fd(i) = (tcm_dv_each1dx(7) - tcm_dv_each(7))/dx;
    dtcm8_fd(i) = (tcm_dv_each1dx(8) - tcm_dv_each(8))/dx;
    dtcmR_fd(i) = (tcmR_dv_total1dx-tcmR_dv_total)/dx;
    dtcmV_fd(i) = (tcm_dv_each1dx(end)-tcm_dv_each(end))/dx;

    dtcm_total_fd(i) = (min_tcm_dv1dx - min_tcm_dv)./dx;

    if i == 15*7+1
        bpp=1;
    end



%     stmCk0_1dx = stm_t1dx(:,:,t1dx==tcm_time(k));
%     stm0Ck_1dx = - J * stmCk0_1dx' * J;
% 
%     stmNCk1dx = stmN01dx * stm0Ck_1dx; % this is including mods to stmN0 from segment 1, when the correction is in seg 2, for example, which is incorrect
%     tCk_idx = find(t1dx==tcm_time(k));

%     stmNCk1dx = dynCellCombine(t1dx, t_s1dx, tCk_idx, target_idx1dx, simparams, stm_t_i1dx);

%     if k == 1
%         tClast_idx = 1;
%     else
%         tClast_idx = find(t==tcm_time(k-1));
%     end

    

%     dstmNCk_fd(:,:,i) = (stmNCk1dx - stmNCk)./dx;


    target_idx1dx = find(t1dx == target_time);        
    tcmk_idx = find(event_times(k) == t1dx);
    stmNCk1dx = dynCellCombine(t1dx, t_s1dx, tcmk_idx, target_idx1dx, simparams, stm_t_i1dx);       
    Tk1dx = [-inv( stmNCk1dx(1:3,4:6) ) * stmNCk1dx(1:3,1:3), -eye(3)];



    Pck_minus1dx = P_i_minus1dx(:,:,k);
%     Pck_minus = P_i_minus(:,:,find(event_times==tcm_time(k)));


    Pc1_minus1dx = P_i_minus1dx(:,:,1);
    Pc2_minus1dx = P_i_minus1dx(:,:,2);
    Pc3_minus1dx = P_i_minus1dx(:,:,3);
    Pc4_minus1dx = P_i_minus1dx(:,:,4);
%     Pc5_minus1dx = P_i_minus1dx(:,:,5);
%     Pcn_minus1dx = P_i_minus1dx(:,:,num_events+1);
%     Pcn_minus1dx = P_i_minus1dx(:,:,num_events);



    dPc1_minusdx_fd(:,:,i) = (Pc1_minus1dx - P_i_minus(:,:,1))./dx;
    dPc2_minusdx_fd(:,:,i) = (Pc2_minus1dx - P_i_minus(:,:,2))./dx;
    dPc3_minusdx_fd(:,:,i) = (Pc3_minus1dx - P_i_minus(:,:,3))./dx;
    dPc4_minusdx_fd(:,:,i) = (Pc4_minus1dx - P_i_minus(:,:,4))./dx;
%     dPc5_minusdx_fd(:,:,i) = (Pcn_minus1dx - P_i_minus(:,:,5))./dx;
%     dPcn_minusdx_fd(:,:,i) = (Pcn_minus1dx - P_i_minus(:,:,num_events))./dx;




    sig_dvk1dx = sqrt(trace(Tk1dx * Pck_minus1dx * Tk1dx'));

    dsig_dvk_dx1dx_fd(i) = (sig_dvk1dx - sig_dvk)./dx;

    dTkdxi_fd(:,:,i) = (Tk1dx - Tk)./dx;

    Ptcm1dx = Tk1dx * Pc1_minus1dx * Tk1dx';

    dPtcm_fd(:,:,i) = (Ptcm1dx - Ptcm)./dx;





    DV1dx = deltaVs_nom1dx(:,nomDVindex);
    DV_mag1dx = vecnorm(DV1dx);
    i_DV1dx = DV1dx/DV_mag1dx;
    
    sig_corrected_DV1dx = i_DV1dx' * Ptcm1dx * i_DV1dx;

    dsig_corrDV_fd = (sig_corrected_DV1dx - sig_corrected_DV)./dx;








    % testing dstmCkClast_dxi


%     if k == 1
%         tClast_idx1dx = 1;
%     else
%         tClast_idx1dx = find(t1dx==tcm_time1dx(k-1));
%     end
%     tCk_idx1dx = find(t1dx==tcm_time1dx(k));
%     stmCkClast1dx = dynCellCombine(t1dx, t_s1dx, tClast_idx1dx, tCk_idx1dx, simparams, stm_t_i1dx);


%     dstmCkClast_dxi(:,:,i) = (stmCkClast1dx - stmCkClast)./dx;


end



% dT_fd_function = T_partial(stmNCk, dstmNCk_fd(:,:,(k-1)*7+1:k*7-1))

% % Get the time the kth correction occurs
% tCk = tcm_time(k);
% 
% % Get the index tCk occurs
% tCk_idx = find(t==tCk);
% 
% % Test which segments have a correction
% tCk_seg = t_s(tCk_idx);
% i=1;
% dstmNCkdxi = calc_dstmNC(x, x_i_f, t, t_s, stm_t, stm_i, stt_t_i, tCk, tCk_idx, tCk_seg, stmN0, i, simparams);
% dTkdxi = T_partial(stmNCk, dstmNCkdxi);
% dsig_dvk_dx1dx_an = calc_dSigkdxi(Tk, dTkdxi,  P_i_minus(:,:,k), dPCkminusdxi);

% [dtcm1_fd' ,  dtcm2_fd' ,  dtcm3_fd', dtcm4_fd', dtcmV_fd]
[dtcm1_fd' ,  dtcm2_fd' ,  dtcm3_fd', dtcmV_fd]


% easier to debug into the calc_multiple_tcm_gradient function and directly
% observe the tcm_gradient_r values

[tcm_gradient, tcm_gradient_r, tcm_gradient_v] = calc_multiple_tcm_gradient(x, x_t, x_i_f, stm_i, stt_i, stm_t, stm_t_i, stt_t_i, t, t_s, tcm_time, tcm_idx, P_i_minus, P_i_plus, deltaVs_nom, simparams);




[repmat((1:7)',simparams.n,1), dtcm1_fd'./tcm_gradient_r(:,1), dtcm1_fd', tcm_gradient_r(:,1)]
[repmat((1:7)',simparams.n,1), dtcm2_fd'./tcm_gradient_r(:,2), dtcm2_fd', tcm_gradient_r(:,2)]
[repmat((1:7)',simparams.n,1), dtcm3_fd'./tcm_gradient_r(:,3), dtcm3_fd', tcm_gradient_r(:,3)]
[repmat((1:7)',simparams.n,1), dtcm4_fd'./tcm_gradient_r(:,4), dtcm4_fd', tcm_gradient_r(:,4)]

[repmat((1:7)',simparams.n,1), dtcm5_fd'./tcm_gradient_r(:,5), dtcm5_fd', tcm_gradient_r(:,5)]
[repmat((1:7)',simparams.n,1), dtcm6_fd'./tcm_gradient_r(:,6), dtcm6_fd', tcm_gradient_r(:,6)]
[repmat((1:7)',simparams.n,1), dtcm7_fd'./tcm_gradient_r(:,7), dtcm7_fd', tcm_gradient_r(:,7)]
[repmat((1:7)',simparams.n,1), dtcm8_fd'./tcm_gradient_r(:,8), dtcm8_fd', tcm_gradient_r(:,8)]



[dtcmV_fd, tcm_gradient_v]
[repmat((1:7)',simparams.n,1), dtcm_total_fd./tcm_gradient, dtcm_total_fd, tcm_gradient]



Mv = [zeros(3,3), eye(3,3)];


for i = 1:simparams.n
    tcm_gradient_v_test((i-1)*7+1:(i-1)*7+6,1) = calc_dSigk(Mv, zeros(3,6,6), P_i_minus(:,:,5), dPcn_minusdx_fd(:,:,(i-1)*7+1:(i-1)*7+6));
    tcm_gradient_v_test((i-1)*7+7,1) = calc_dSigk(Mv, zeros(3,6), P_i_minus(:,:,5), dPcn_minusdx_fd(:,:,(i-1)*7+7));
%     calc_dSigk(Mv, zeros(3,6,6), P_i_minus(:,:,5), dPcn_minusdx_fd(:,:,162:167))
%     calc_dSigk(Mv, zeros(3,6), P_i_minus(:,:,5), dPcn_minusdx_fd(:,:,168))

end



% when i=2 (second segment) and k=2 (second correction), verifying dTdxi and
% dstmCkClastdxi (stmC2C1)


