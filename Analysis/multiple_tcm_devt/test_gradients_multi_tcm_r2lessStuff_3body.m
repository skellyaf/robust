%% Multiple TCM stuff
% reference derivation in lyx document: stm_target_vel_correction

clear;
clear global;
close all; 
clc;
format longg;
cd('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust_2body_tcmerror');
addpath(genpath('./'));

%%
% init_fn = 'C:\Users\skell\OneDrive - USU\Documents\code_repos\robust_2body_tcmerror\init_traj_files\initialize_simulation_parameters_leo28_to_geo0';
init_fn = './init_traj_files/init_simparams_cr3bp_leo_to_llo';



run(init_fn);
% simparams = generateInitialXfer(simparams);
%% ONLY If using a previous reference trajectory as the initial guess:
load('.\init_traj_files\initial_guesses\tli_llo_deterministic_opt.mat')
% load('.\init_traj_files\initial_guesses\polar_llo_to_nrho_apolune.mat')
simparams.x0 = x_opt;

%%

x = simparams.x0;
dx = -sqrt(eps);

% Symplectic unit matrix
J = [zeros(3,3), eye(3,3); -eye(3,3), zeros(3,3)];

[stm_i, stt_i, x_i_f, x_t, stm_t, stt_t_i, t, t_s, stm_t_i]  = createStateStmSttHistory(x, simparams);

% % Using saved dynamics, calculate total impulsive delta V
% [deltaV, deltaVs_nom] = calcDeltaV(x, x_i_f, simparams);

% Using saved dynamics, calculate the TCMs

[tcm_time,tcm_idx,min_tcm_dv] = opt_multiple_tcm(x, t, stm_t, simparams);

% 
[~, ~, tcm_dv_each, P_i_minus, P_i_plus] = calc_covariance_tcmdv(x, t, stm_t, tcm_time, simparams);

%%
R = simparams.R;
G = [zeros(3,3); eye(3,3)];

x = reshape(x,simparams.m,simparams.n);


if simparams.target_final_maneuver
    final_seg = simparams.maneuverSegments(end);
    target_time = sum(x(7,1:final_seg - 1));
    target_idx = find(target_time == t) ;
    stmN0 = stm_t(:,:,target_idx);
    stm_t = stm_t(:,:,1:target_idx);
else
    stmN0 = stm_t(:,:,end);
end



k=4;
% stmCk0 = stm_t(:,:,tcm_idx(k));
% stm0Ck = -J * stmCk0' * J;
% stmNCk = stmN0 * stm0Ck;
% target_idx=360;
stmNCk = dynCellCombine(t, t_s, tcm_idx(k), target_idx, simparams, stm_t_i);


Tk = [-inv( stmNCk(1:3,4:6) ) * stmNCk(1:3,1:3), -eye(3)];
Pck_minus = P_i_minus(:,:,k);
sig_dvk = sqrt(trace(Tk * Pck_minus * Tk'));

[~, tcm_dv_total, tcm_dv] = calc_covariance_tcmdv(x, t, stm_t, tcm_time, simparams);
tcmR_dv_total = sum(tcm_dv(1:end-1));
tcm1_dv = tcm_dv(1);
tcm2_dv = tcm_dv(2);
tcm3_dv = tcm_dv(3);
tcm4_dv = tcm_dv(4);
dtcmR_fd = zeros(length(x(:)),1);
dtcmV_fd = zeros(length(x(:)),1);
dtcm_total_fd=zeros(length(x(:)),1);
dTkdxi = zeros(3,6,length(x(:)));
dstmNCk_fd = zeros(6,6,length(x(:)));


t20_idx = find(t==x(7,1));
t2f_idx = find(t==sum(x(7,1:2)));
[stm_20_tc2, stt_20_tc2]  = dynCellCombine(t, t_s, tcm_idx(k), t20_idx, simparams, stm_t_i, stt_t_i);
stm_2f_20  = dynCellCombine(t, t_s, t20_idx, t2f_idx, simparams, stm_t_i);



if k == 1
    tClast_idx = 1;
else
    tClast_idx = find(t==tcm_time(k-1));
end
tCk_idx = find(t==tcm_time(k));
stmCkClast = dynCellCombine(t, t_s, tClast_idx, tCk_idx, simparams, stm_t_i);



parfor i = 1:simparams.m*simparams.n
    i
    x1dx = x(:);
    x1dx(i) = x(i) + dx;
    [stm_i1dx, stt_i1dx, x_i_f1dx, x_t1dx, stm_t1dx, stt_t_i1dx, t1dx, t_s1dx, stm_t_i1dx]  = createStateStmSttHistory(x1dx, simparams);
%     [x_t1dx, stm_t1dx, stt_t_i1dx, t1dx, t_s1dx] = addToStateStmSttHistory(x_t1dx, stm_t1dx, stt_t_i1dx, t1dx, t_s1dx, [tcm_time, target_time], simparams);
    [x_t1dx, stm_t1dx, stm_t_i1dx stt_t_i1dx, t1dx, t_s1dx] = addToStateStmSttHistory(x_t1dx, stm_t1dx, stm_t_i1dx, stt_t_i1dx, t1dx, t_s1dx, [tcm_time], simparams);

    if simparams.target_final_maneuver
        final_seg = simparams.maneuverSegments(end);
        x1dx = reshape(x1dx,simparams.m,simparams.n);
        target_time = sum(x1dx(7,1:final_seg - 1));
        target_idx1dx = find(target_time == t1dx) ;
        stmN01dx = stm_t1dx(:,:,target_idx1dx);
        stm_t1dx = stm_t1dx(:,:,1:target_idx1dx);
        t1dx = t1dx(1:target_idx1dx);
    else
        stmN01dx = stm_t1dx(:,:,end);
    end

%     if x1dx(7) + x1dx(14) ~= target_time
%         x1dx(14) = target_time - x1dx(7);
%     end

    [tcm_time1dx,tcm_idx1dx,min_tcm_dv1dx,~,~,tcm_num_option_DVs1dx] = opt_multiple_tcm(x1dx, t1dx, stm_t1dx, simparams);
%     [tcm_time,tcm_idx,min_tcm_dv,tcm_time_cell,tcm_idx_cell,tcm_num_option_DVs] 
    [~, tcm_dv_total1dx, tcm_dv1dx] = calc_covariance_tcmdv(x1dx, t1dx, stm_t1dx, tcm_time, simparams);
    [~, ~, tcm_dv_each1dx, P_i_minus1dx, P_i_plus1dx] = calc_covariance_tcmdv(x1dx, t1dx, stm_t1dx, tcm_time, simparams);

    tcmR_dv_total1dx = sum(tcm_dv1dx(1:end-1));

    dtcm1_fd(i) = (tcm_dv1dx(1) - tcm_dv(1))/dx;
    dtcm2_fd(i) = (tcm_dv1dx(2) - tcm_dv(2))/dx;
    dtcm3_fd(i) = (tcm_dv1dx(3) - tcm_dv(3))/dx;
    dtcm4_fd(i) = (tcm_dv1dx(4) - tcm_dv(4))/dx;
    dtcmR_fd(i) = (tcmR_dv_total1dx-tcmR_dv_total)/dx;
    dtcmV_fd(i) = (tcm_dv1dx(end)-tcm_dv(end))/dx;

    dtcm_total_fd(i) = (tcm_dv_total1dx - tcm_dv_total)./dx;

    



%     stmCk0_1dx = stm_t1dx(:,:,t1dx==tcm_time(k));
%     stm0Ck_1dx = - J * stmCk0_1dx' * J;
% 
%     stmNCk1dx = stmN01dx * stm0Ck_1dx; % this is including mods to stmN0 from segment 1, when the correction is in seg 2, for example, which is incorrect
    tCk_idx = find(t1dx==tcm_time(k));
    stmNCk1dx = dynCellCombine(t1dx, t_s1dx, tCk_idx, target_idx1dx, simparams, stm_t_i1dx);

    if k == 1
        tClast_idx = 1;
    else
        tClast_idx = find(t==tcm_time(k-1));
    end

    

    dstmNCk_fd(:,:,i) = (stmNCk1dx - stmNCk)./dx;

    Tk1dx = [-inv( stmNCk1dx(1:3,4:6) ) * stmNCk1dx(1:3,1:3), -eye(3)];

    Pck_minus1dx = P_i_minus1dx(:,:,k);

    Pc1_minus1dx = P_i_minus1dx(:,:,1);
    Pc2_minus1dx = P_i_minus1dx(:,:,2);
    Pc3_minus1dx = P_i_minus1dx(:,:,3);
    Pc4_minus1dx = P_i_minus1dx(:,:,4);
    Pcn_minus1dx = P_i_minus1dx(:,:,5);

    dPc1_minusdx_fd(:,:,i) = (Pc1_minus1dx - P_i_minus(:,:,1))./dx;
    dPc2_minusdx_fd(:,:,i) = (Pc2_minus1dx - P_i_minus(:,:,2))./dx;
    dPc3_minusdx_fd(:,:,i) = (Pc3_minus1dx - P_i_minus(:,:,3))./dx;
    dPc4_minusdx_fd(:,:,i) = (Pc4_minus1dx - P_i_minus(:,:,4))./dx;
    dPcn_minusdx_fd(:,:,i) = (Pcn_minus1dx - P_i_minus(:,:,5))./dx;




    sig_dvk1dx = sqrt(trace(Tk1dx * Pck_minus1dx * Tk1dx'));

    dsig_dvk_dx1dx_fd(i) = (sig_dvk1dx - sig_dvk)./dx;

    dTkdxi(:,:,i) = (Tk1dx - Tk)./dx;



    % testing dstmCkClast_dxi


    if k == 1
        tClast_idx1dx = 1;
    else
        tClast_idx1dx = find(t1dx==tcm_time1dx(k-1));
    end
    tCk_idx1dx = find(t1dx==tcm_time1dx(k));
    stmCkClast1dx = dynCellCombine(t1dx, t_s1dx, tClast_idx1dx, tCk_idx1dx, simparams, stm_t_i1dx);


    dstmCkClast_dxi(:,:,i) = (stmCkClast1dx - stmCkClast)./dx;


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

[dtcm1_fd' ,  dtcm2_fd' ,  dtcm3_fd', dtcm4_fd', dtcmV_fd]


% easier to debug into the calc_multiple_tcm_gradient function and directly
% observe the tcm_gradient_r values
tic
[tcm_gradient, tcm_gradient_r, tcm_gradient_v] = calc_multiple_tcm_gradient(x, x_t, x_i_f, stm_i, stt_i, stm_t, stm_t_i, stt_t_i, t, t_s, tcm_time, tcm_idx, tcm_dv_each, P_i_minus, P_i_plus, simparams);
toc

[dtcmV_fd, tcm_gradient_v]
[dtcm_total_fd, tcm_gradient]



Mv = [zeros(3,3), eye(3,3)];


for i = 1:simparams.n
    tcm_gradient_v_test((i-1)*7+1:(i-1)*7+6,1) = calc_dSigk(Mv, zeros(3,6,6), P_i_minus(:,:,5), dPcn_minusdx_fd(:,:,(i-1)*7+1:(i-1)*7+6));
    tcm_gradient_v_test((i-1)*7+7,1) = calc_dSigk(Mv, zeros(3,6), P_i_minus(:,:,5), dPcn_minusdx_fd(:,:,(i-1)*7+7));
%     calc_dSigk(Mv, zeros(3,6,6), P_i_minus(:,:,5), dPcn_minusdx_fd(:,:,162:167))
%     calc_dSigk(Mv, zeros(3,6), P_i_minus(:,:,5), dPcn_minusdx_fd(:,:,168))

end



% when i=2 (second segment) and k=2 (second correction), verifying dTdxi and
% dstmCkClastdxi (stmC2C1)


