%% Multiple TCM stuff
% reference derivation in lyx document: stm_target_vel_correction

clear;
clear global;
close all; 
clc;
format longg;
addpath(genpath('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust'));
% addpath(genpath('../'));
% addpath('C:\Users\skell\OneDrive - USU\Documents\PhD\research\directional state transition tensors');

savename = ['testing_dv_stats_nri'];
saveOutput = true; % bool for saving the output or not, true or false
saveVideo = false;


%% Create simulation parameters structure by running initialization script

% init_fn = './init_traj_files/init_simparams_cr3bp_leo_lloflyby_nri_3dv';
% init_fn = 'init_simparams_cr3bp_leoinclined_lloflyby_nri_3dv';
init_fn = './init_traj_files/init_simparams_cr3bp_nrho_rdvz_2dv';

run(init_fn);

%% ONLY If using a previous reference trajectory as the initial guess:
% load('.\init_traj_files\initial_guesses\tli_llo_deterministic_opt.mat')
% load('.\init_traj_files\initial_guesses\polar_llo_to_nrho_apolune.mat')

% load('eed_leo_planar_13day.mat');
% load('nri_det_opt.mat');
% load('nri_planar_det_opt.mat');
% load('leo_plf_dro3_detOpt.mat')
load('nrho_apolune_rdvz.mat');

simparams.x0 = x_opt;


%% Initialize transfer segments
% simparams = generateInitialXfer(simparams); % only for 2bp

%% View problem setup


% [~, x_i_f0, x_t0, stm_t0, t0, t_s0] = createStateStmHistory(simparams.x0, simparams);
% [stm_i0, stt_i0, x_i_f0, x_t0, stm_t0, stt_t_i, t0, t_s0, stm_t_i0]  = createStateStmSttHistory(simparams.x0, simparams);
[traj]  = createStateStmSttQdQHistory(simparams.x0, simparams);

% tcm_idx = [100, 1000, 2000, 5500, 6000];
% tcm_time = traj.t(tcm_idx)';
[deltaV, deltaVs_nom] = calcDeltaV(simparams.x0, traj.x_i_f, traj.stm_i, simparams);

[tcm_time, tcm_idx, min_tcm_dv, ~, ~, tcm_dv_each] = opt_multiple_tcm_wQ(simparams.x0, traj, deltaVs_nom, simparams); % inputs: x, t, t_s, stm_t, stm_t_i, simparams

figure
plotMultiSegTraj(x_opt, traj.x_t, traj.t_s, simparams, tcm_idx)

Q_k_minus = calc_Q_events(traj, simparams.x0, tcm_time, simparams);


% Calculate total impulsive delta V for initial guess trajectory
% [deltaV0, deltaVs_nom0] = calcDeltaV(simparams.x0,x_i_f0,simparams);
% [tcm_min, tcm_time, tcm_r, tcm_v, ~, ~, ~, ~, tcm_total_t] = tcmPair_rv(x_opt, t, stm_t, deltaVs_nom, simparams);

tic
% [tcm_time0, tcm_idx0, min_tcm_dv0, P_i_minus, P_i_plus, tcm_dv_each0] = opt_multiple_tcm(simparams.x0, deltaVs_nom0, t0, t_s0, stm_t0, stm_t_i0, simparams); % inputs: x, t, t_s, stm_t, stm_t_i, simparams
toc






%%
R = simparams.R;
G = [zeros(3,3); eye(3,3)];
x = simparams.x0;
% x = reshape(x,simparams.m,simparams.n);




%% Also testing calc_covariance_wQ functions here

[event_times, event_indicator] = define_events_v2(x(:), traj.t, tcm_time, simparams);
event_idx_logical = logical(sum(traj.t'==event_times', 1));    
event_idxs = find(event_idx_logical);

[Q_k_km1, dQ_k_km1_dxi, dQ_k_km1_ddti] = calc_Q_events(traj, x, tcm_time, simparams);





% [P, tcm_dv_total, tcm_dv, P_i_minus, P_i_plus] = calc_covariance_wQ_tcmdv(x, traj, tcm_time, 1, deltaVs_nom, P_i, simparams);


[P, tcm_dv_total, tcm_dv, P_i_minus, P_i_plus] = calc_covariance_wQ_tcmdv_v3(x, traj, tcm_time, 1, deltaVs_nom, simparams.P_initial, Q_k_km1, simparams);
% [Pnoq, tcm_dv_totalnoq, tcm_dvnoq, P_i_minusnoq, P_i_plusnoq] = calc_covariance_wQ_tcmdv_v3(x, traj, tcm_time, 1, deltaVs_nom, simparams.P_initial, zeros(6,6,8), simparams);



%% Testing the gradient of Qbar at each event


x=reshape(x,simparams.m,simparams.n);


% dx = sqrt(eps);
dx = 1e-5;
% dx = 1e-10;

dQ_1_minus_dx = zeros(6,6,length(x(:)));
dQ_2_minus_dx = zeros(6,6,length(x(:)));
dQ_3_minus_dx = zeros(6,6,length(x(:)));
dQ_4_minus_dx = zeros(6,6,length(x(:)));
dQ_5_minus_dx = zeros(6,6,length(x(:)));
dQ_6_minus_dx = zeros(6,6,length(x(:)));
dQ_7_minus_dx = zeros(6,6,length(x(:)));
dQ_8_minus_dx = zeros(6,6,length(x(:)));


dtcm_1_fd = zeros(length(x(:)),1);
dtcm_2_fd = zeros(length(x(:)),1);
dtcm_3_fd = zeros(length(x(:)),1);
% dtcm_4_fd = zeros(length(x(:)),1);
% dtcm_5_fd = zeros(length(x(:)),1);
% dtcm_6_fd = zeros(length(x(:)),1);
% dtcm_7_fd = zeros(length(x(:)),1);
% dtcm_8_fd = zeros(length(x(:)),1);

dtcm_total_fd = zeros(length(x(:)),1);
M = [eye(3), zeros(3,3)];
Pr_constraint = simparams.P_max_r - sqrt(trace(M*P*M'));


parfor i = 1:simparams.m*simparams.n
% for i = 1:simparams.m*simparams.n
    i
    x1dx = x(:);
    x1dx(i) = x(i) + dx;
    modseg = floor((i-1)/7)+1;
    [traj1dx]  = createStateStmSttQdQHistory(x1dx, simparams);
%     [traj1dx.x_t, traj1dx.stm_t, traj1dx.stm_t_i, traj1dx.stt_t_i, traj1dx.t, traj1dx.t_s, traj1dx.Q_t, traj1dx.Q_t_i] = addToStateStmSttQHistory(traj1dx.x_t, traj1dx.stm_t, traj1dx.stm_t_i, traj1dx.stt_t_i, traj1dx.t, traj1dx.t_s, traj1dx.Q_t, traj1dx.Q_t_i, [event_times], simparams);





%%%%% tcm time testing

%     if i == 28
%         ppp=1;  
%     end
% 
% 
%     % calculate time after the end of segment 5 / beginning of segment 6
%     % that the TCM occurs (pre FD)
%     k = 2; % the 2nd TCM
% 
%     time_k = tcm_time(k);
%     tcm_seg = traj.t_s(tcm_idx(k));
% 
%     x1dx = reshape(x1dx,simparams.m,simparams.n);
%     endModSegTotalTime = sum(x(7,1:modseg));
%     endModSegTotalTime1dx = sum(x1dx(7,1:modseg));
% 
% 
%     idxs_tcmseg = find(traj.t_s == tcm_seg);
%     dt_past_modseg = time_k - endModSegTotalTime;
% 
% 
% 
%     tcm_idx1dx = find(tcm_time(k) == traj1dx.t)
%     dt_past_modseg1dx = traj1dx.t(tcm_idx1dx) - endModSegTotalTime1dx
%     dt_past_modseg1dx_plus1 = traj1dx.t(tcm_idx1dx+1) - endModSegTotalTime1dx







%%%%%%%%%%%%%


%     x1dx = x1dx(:);






    if mod(i,7) == 0
        


        tcm_seg = traj.t_s(tcm_idx);
%         tcm_time1dx = tcm_time;
% 
%         for k = find(tcm_seg > modseg,1) : length(tcm_seg)
%             kidx = find(tcm_time(k) + dx == traj1dx.t);
%             
%             tcm_time1dx(k) = traj1dx.t(kidx);
% 
%             assert(tcm_time1dx(k) == tcm_time(k) + dx);
% 
%         end


        modtcmidx = find(tcm_seg > modseg);
        tcm_time1dx = tcm_time;
        tcm_time1dx(modtcmidx) = tcm_time(modtcmidx) + dx;

    else
        tcm_time1dx = tcm_time;
    end


    
    

    [traj1dx.x_t, traj1dx.stm_t, traj1dx.stm_t_i, traj1dx.stt_t_i, traj1dx.t, traj1dx.t_s, traj1dx.Q_t, traj1dx.Q_t_i] = addToStateStmSttQHistory(traj1dx.x_t, traj1dx.stm_t, traj1dx.stm_t_i, traj1dx.stt_t_i, traj1dx.t, traj1dx.t_s, traj1dx.Q_t, traj1dx.Q_t_i, [event_times, tcm_time1dx], simparams);






    Q_k_minus1dx = calc_Q_events(traj1dx, x1dx, tcm_time1dx, simparams);

    [deltaV1dx, deltaVs_nom1dx] = calcDeltaV(x1dx, traj1dx.x_i_f, traj1dx.stm_i, simparams);


    [P1dx, tcm_dv_total1dx, tcm_dv1dx, P_i_minus1dx, P_i_plus1dx] = calc_covariance_wQ_tcmdv_v3(x1dx, traj1dx, tcm_time1dx, 1, deltaVs_nom1dx, simparams.P_initial, Q_k_minus1dx, simparams);



    dQ_1_minus_dx(:,:,i) = (Q_k_minus1dx(:,:,1) - Q_k_minus(:,:,1)) ./ dx;
    dQ_2_minus_dx(:,:,i) = (Q_k_minus1dx(:,:,2) - Q_k_minus(:,:,2)) ./ dx; 
    dQ_3_minus_dx(:,:,i) = (Q_k_minus1dx(:,:,3) - Q_k_minus(:,:,3)) ./ dx;
    dQ_4_minus_dx(:,:,i) = (Q_k_minus1dx(:,:,4) - Q_k_minus(:,:,4)) ./ dx;
    dQ_5_minus_dx(:,:,i) = (Q_k_minus1dx(:,:,5) - Q_k_minus(:,:,5)) ./ dx;
%     dQ_6_minus_dx(:,:,i) = (Q_k_minus1dx(:,:,6) - Q_k_minus(:,:,6)) ./ dx;
%     dQ_7_minus_dx(:,:,i) = (Q_k_minus1dx(:,:,7) - Q_k_minus(:,:,7)) ./ dx;
%     dQ_8_minus_dx(:,:,i) = (Q_k_minus1dx(:,:,8) - Q_k_minus(:,:,8)) ./ dx;


    % Calculations

    dtcm_1_fd(i) = (tcm_dv1dx(1) - tcm_dv(1)) ./ dx;
    dtcm_2_fd(i) = (tcm_dv1dx(2) - tcm_dv(2)) ./ dx;
    dtcm_3_fd(i) = (tcm_dv1dx(3) - tcm_dv(3)) ./ dx;
    dtcm_4_fd(i) = (tcm_dv1dx(4) - tcm_dv(4)) ./ dx;
%     dtcm_5_fd(i) = (tcm_dv1dx(5) - tcm_dv(5)) ./ dx;
%     dtcm_6_fd(i) = (tcm_dv1dx(6) - tcm_dv(6)) ./ dx;
%     dtcm_7_fd(i) = (tcm_dv1dx(7) - tcm_dv(7)) ./ dx;
%     dtcm_8_fd(i) = (tcm_dv1dx(8) - tcm_dv(8)) ./ dx;

    dtcm_total_fd(i) = (tcm_dv_total1dx - tcm_dv_total) ./ dx;


    

        % Testing covariance matrices
    P1_fd(:,:,i) = (P_i_minus1dx(:,:,1) - P_i_minus(:,:,1)) ./ dx;
    P2_fd(:,:,i) = (P_i_minus1dx(:,:,2) - P_i_minus(:,:,2)) ./ dx;
    P3_fd(:,:,i) = (P_i_minus1dx(:,:,3) - P_i_minus(:,:,3)) ./ dx;
    P4_fd(:,:,i) = (P_i_minus1dx(:,:,4) - P_i_minus(:,:,4)) ./ dx;
    P5_fd(:,:,i) = (P_i_minus1dx(:,:,5) - P_i_minus(:,:,5)) ./ dx;
%     P6_fd(:,:,i) = (P_i_minus1dx(:,:,6) - P_i_minus(:,:,6)) ./ dx;
%     P7_fd(:,:,i) = (P_i_minus1dx(:,:,7) - P_i_minus(:,:,7)) ./ dx;


    Pr_constraint1dx = simparams.P_max_r-sqrt(trace(M*P1dx*M'));

    dPr_constraint(i) = (Pr_constraint1dx - Pr_constraint) / dx;


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


% save('Qbar_gradients.mat')
% save('tcm_wQ_gradients.mat')

%% Gradient comparison

t_idx = 7:7:simparams.m*simparams.n;



[P, tcm_dv_total, tcm_dv, P_i_minus, P_i_plus] = calc_covariance_wQ_tcmdv_v3(x, traj, tcm_time, 1, deltaVs_nom, simparams.P_initial, Q_k_km1, simparams);

[Q_k_km1, dQ_k_km1_dxi, dQ_k_km1_ddti] = calc_Q_events(traj, x, tcm_time, simparams);


[tcm_gradient, tcm_gradient_r, tcm_gradient_v, dPCkminusdxi, dPCkminusddti] = calc_multiple_tcm_gradient_wQ(x, traj, tcm_time, tcm_idx, P_i_minus, dQ_k_km1_dxi, dQ_k_km1_ddti, deltaVs_nom, simparams);

%%%% POTENTIAL BUG
% Do some debugging: The analytical dP for the 5th TCM event (the target)
% minus, is coming up as zeros. Don't think that should be the case because
% that is the target covariance (minus) and modifications to the trajectory
% should be impacting the target covariance. Tried modifying the start P
% growth node, maybe there was a misnumbering bug...no difference. 

% Comparing the P fd and P analytical

% dP dt

P5_fd(:,:,t_idx)
dPCkminusddti(:,:,5,:) % is making all zeros
dPCkminusddti(:,:,4,:) % is not


% These aren't matching
P5_fd(:,:,t_idx(7))
dPCkminusddti(:,:,5,7)

% The following are matching
P4_fd(:,:,t_idx(7))
dPCkminusddti(:,:,4,7)









% Analytical Pr constraint gradient
for i = 1:simparams.n
    for j = 1:6
        dPr_an((i-1)*7 + j) = 1/2 * trace(M*P*M')^(-1/2) * -trace(M*dPCkminusdxi(:,:,j,5,i)*M');
    end
        dPr_an(i*7) = 1/2 * trace(M*P*M')^(-1/2) * -trace(M*dPCkminusddti(:,:,5,i)*M');
end

[dPr_constraint', [dPr_an]']



% t0i = tkm1, tk < tfi
dQ_1_minus_dx(:,:,t_idx)
dQ_k_km1_ddti(:,:,1,:)

% t0_i < tkm1, tk = tfi (% CASE 1)
dQ_2_minus_dx(:,:,t_idx) 
dQ_k_km1_ddti(:,:,2,:) % matches


dQ_3_minus_dx(:,:,t_idx) % issue with this set...first 4 matrices both exist, but don't match 
dQ_k_km1_ddti(:,:,3,:)

dQ_4_minus_dx(:,:,t_idx) % has matrices 1-9 populated
dQ_k_km1_ddti(:,:,4,:) % only has matrices 5-9 populated...and numbers aren't matching

dQ_5_minus_dx(:,:,t_idx) % has 1-13
dQ_k_km1_ddti(:,:,5,:) % only has 10-13

dQ_6_minus_dx(:,:,t_idx) % has 1-14
dQ_k_km1_ddti(:,:,6,:) % % has 13-19

dQ_7_minus_dx(:,:,t_idx) % none populated
dQ_k_km1_ddti(:,:,7,:) % 20 is populated

dQ_8_minus_dx(:,:,t_idx)
dQ_k_km1_ddti(:,:,8,:) % 21-24 is populated




[dtcm_1_fd, tcm_gradient_r(:,1)]
[dtcm_2_fd, tcm_gradient_r(:,2), dtcm_total_fd-tcm_gradient]
[dtcm_3_fd, tcm_gradient_r(:,3)]
[dtcm_4_fd, tcm_gradient_r(:,4), dtcm_total_fd-tcm_gradient]
[dtcm_5_fd, tcm_gradient_r(:,5), dtcm_total_fd-tcm_gradient]
[dtcm_6_fd, tcm_gradient_r(:,6), dtcm_total_fd-tcm_gradient]
[dtcm_7_fd, tcm_gradient_r(:,7), dtcm_total_fd-tcm_gradient]
[dtcm_8_fd, tcm_gradient_v(:,1)]

[dtcm_total_fd, tcm_gradient, dtcm_total_fd-tcm_gradient]
[dtcm_total_fd, tcm_gradient, (dtcm_total_fd-tcm_gradient)./dtcm_total_fd]





%% Gradient development
% load('Qbar_gradients.mat')
% load('tcm_wQ_gradients.mat')

x=reshape(x,simparams.m,simparams.n);
k = 2;
tcm_k_time = tcm_time(k);
tcm_k_idx = tcm_idx(k);
tcm_k_seg = traj.t_s(tcm_k_idx);


if k > 1
    km1 = k - 1;
    tcm_km1_time = tcm_time(km1);
    tcm_km1_idx = tcm_idx(km1);
    tcm_km1_seg = traj.t_s(tcm_km1_idx);
else
    tcm_km1_time = 0;
    tcm_km1_idx = 1;
    tcm_km1_seg = 1;
end

if tcm_km1_seg == 1
    ti = 0;
    i_idx = 1;
    i_seg = 1;
else
    ti = sum(x(7,1:tcm_km1_seg-1));
    i_idx = find(traj.t == ti);
    i_seg = t_s(i_idx);
end



t_fi = sum(x(7,1:tcm_km1_seg));
if_idx = find(traj.t == t_fi);



    
tic
[Q_k_km1, dQ_k_km1_dxi] = calc_Q_events(traj, x, tcm_time, simparams);
toc


k=7
i_seg=15
dQ_k_km1_dxi(:,:,:,k,i_seg)


%%%% comparing
dQ_1_minus_dx(:,:,i_seg*7-6:i_seg*7-1) 
dQ_2_minus_dx(:,:,i_seg*7-6:i_seg*7-1) 
dQ_3_minus_dx(:,:,i_seg*7-6:i_seg*7-1)
dQ_4_minus_dx(:,:,i_seg*7-6:i_seg*7-1) 
dQ_5_minus_dx(:,:,i_seg*7-6:i_seg*7-1)
dQ_6_minus_dx(:,:,i_seg*7-6:i_seg*7-1)
dQ_7_minus_dx(:,:,i_seg*7-6:i_seg*7-1)




%%%%%%%%%% COME BACK REMINDER: FIX THE OPT TCM ALGORITHM - HAVE BEEN HAVING
%%%%%%%%%% ISSUE WITH UNINTENTIONAL TCM & NOM MANEUVERS BEING TESTING ON
%%%%%%%%%% TOP OF EACH OTHER


%% leftover


[dtcm1_fd' ,  dtcm2_fd' ,  dtcm3_fd', dtcmV_fd]


% easier to debug into the calc_multiple_tcm_gradient function and directly
% observe the tcm_gradient_r values

[tcm_gradient, tcm_gradient_r, tcm_gradient_v] = calc_multiple_tcm_gradient_wQ(x, traj.x_t, traj.x_i_f, traj.stm_i, traj.stt_i, traj.stm_t, traj.stm_t_i, traj.stt_t_i, traj.t, traj.t_s, tcm_time, tcm_idx, P_i_minus, dQ_k_km1_dxi, deltaVs_nom, simparams);
[tcm_gradient, tcm_gradient_r, tcm_gradient_v] = calc_multiple_tcm_gradient_wQ(x, traj.x_t, traj.x_i_f, traj.stm_i, traj.stt_i, traj.stm_t, traj.stm_t_i, traj.stt_t_i, traj.t, traj.t_s, tcm_time, tcm_idx, P_i_minus, zeros(6,6,6,8,simparams.m*simparams.n), deltaVs_nom, simparams);

[dtcm_1_fd, tcm_gradient_r(:,1)]
[dtcm_2_fd, tcm_gradient_r(:,2)]
[dtcm_3_fd, tcm_gradient_r(:,3)]
[dtcm_4_fd, tcm_gradient_r(:,4)]




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


