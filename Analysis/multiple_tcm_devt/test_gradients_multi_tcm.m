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
init_fn = 'C:\Users\skell\OneDrive - USU\Documents\code_repos\robust_2body_tcmerror\init_traj_files\initialize_simulation_parameters_leo28_to_geo0';
run(init_fn);
simparams = generateInitialXfer(simparams);

%%

x = simparams.x0;

[stm_i, stt_i, x_i_f, x_t, stm_t, stt_t_i, t, t_s, stm_t_i]  = createStateStmSttHistory(x, simparams);

% % Using saved dynamics, calculate total impulsive delta V
% [deltaV, deltaVs_nom] = calcDeltaV(x, x_i_f, simparams);

% Using saved dynamics, calculate the TCMs

[tcm_time,tcm_idx,min_tcm_dv] = opt_multiple_tcm(x, t, stm_t, simparams);

% 
[~, ~, tcm_dv_each, P_i_minus, P_i_plus] = calc_covariance_tcmdv(x, t, stm_t, tcm_time, simparams);

%% Finite difference calculation of P_1_minus_dx1

dx = sqrt(eps);
dPc1minusdx1_fd = zeros(6,6,6);


for i = 1:6
    x1dx = x;
    x1dx(i) = x(i) + dx;
    [stm_i1dx, stt_i1dx, x_i_f1dx, x_t1dx, stm_t1dx, stt_t_i1dx, t1dx, t_s1dx]  = createStateStmSttHistory(x1dx, simparams);
    [tcm_time1dx,tcm_idx1dx,min_tcm_dv1dx] = opt_multiple_tcm(x1dx, t1dx, stm_t1dx, simparams);
    [~, ~, tcm_dv_each1dx, P_i_minus1dx, P_i_plus1dx] = calc_covariance_tcmdv(x1dx, t1dx, stm_t1dx, tcm_time1dx, simparams);

    dPc1minusdx1_fd(:,:,i) = (P_i_minus1dx(:,:,1) - P_i_minus(:,:,1))./dx;


end
J = [zeros(3,3), eye(3,3); -eye(3,3), zeros(3,3)];


% Analytical comparison for P_1_minus_dx1
dPc1minusdx1_an = zeros(6,6,6);
stmC10 = stm_t(:,:,tcm_idx(1));
sttC10 = stt_t_i{1}(:,:,:,tcm_idx(1));

P0 = simparams.P_initial;


dPc1minusdx1_an = tmult(stmC10,tmult(P0,sttC10,[0 1])) + tmult(sttC10,tmult(P0,stmC10,[0 1]));
%%%%% works!!!!
% in the function:
dPc1minusdx1_anf = calc_dPckMinus(stmC10, sttC10, 1, simparams);

%% Finite difference for P_2_minus_dx1
dPc2minusdx1_fd = zeros(6,6,21);



x = reshape(x,simparams.m,simparams.n);
% STM from beginning of trajectory to the state being targeted
if simparams.target_final_maneuver
    final_seg = simparams.maneuverSegments(end);
    target_time = sum(x(7,1:final_seg - 1));
    target_idx = find(target_time == t) ;
    stmN0 = stm_t(:,:,target_idx);
    stm_t = stm_t(:,:,1:target_idx);
else
    stmN0 = stm_t(:,:,end);
end


stmC20 = stm_t(:,:,tcm_idx(2));

stm0C1 = -J * stmC10' * J;
stmC2C1 = stmC20 * stm0C1;

dstmC2C1dx1_fd = zeros(6,6,21);



% 
% 
% for i = 1:21
%     x1dx = x(:);
%     x1dx(i) = x(i) + dx;
%     [stm_i1dx, stt_i1dx, x_i_f1dx, x_t1dx, stm_t1dx, stt_t_i1dx, t1dx, t_s1dx]  = createStateStmSttHistory(x1dx, simparams);
%     [tcm_time1dx,tcm_idx1dx,min_tcm_dv1dx] = opt_multiple_tcm(x1dx, t1dx, stm_t1dx, simparams);
%     [x_t1dx, stm_t1dx, stt_t_i1dx, t1dx, t_s1dx] = addToStateStmSttHistory(x_t1dx, stm_t1dx, stt_t_i1dx, t1dx, t_s1dx, [tcm_time], simparams);
%     [~, ~, tcm_dv_each1dx, P_i_minus1dx, P_i_plus1dx] = calc_covariance_tcmdv(x1dx, t1dx, stm_t1dx, tcm_time, simparams);
% 
%     dPc2minusdx1_fd(:,:,i) = (P_i_minus1dx(:,:,2) - P_i_minus(:,:,2))./dx;
% 
% 
% 
%     stmC10_1dx = stm_t1dx(:,:,tcm_idx1dx(1));
%     stm0C1_1dx = - J * stmC10_1dx' * J;
%     stmC20_1dx = stm_t1dx(:,:,tcm_idx1dx(2));
%     stmC2C1_1dx = stmC20_1dx * stm0C1_1dx;
%     dstmC2C1dx1_fd(:,:,i) = (stmC2C1_1dx - stmC2C1)./dx;
% 
% end

% Analytical comparison
dPc2minusdx1_an = zeros(6,6,6);





% Getting sttC2C1
% options: sttC20 * stt0C1
%       sttC2S2 * sttS20 * stt0C1

seg2_c2_idx = tcm_idx(2) - size(stt_t_i{1},4) + 1; %%%%% THIS COULD BE ONE OFF
sttC2S2 = stt_t_i{t_s(tcm_idx(2))}(:,:,:,seg2_c2_idx);

seg2_indices = find(t_s==2);
seg2_start_idx = seg2_indices(1);
stmS20 = stm_t(:,:,seg2_start_idx);
stm0S2 = -J * stmS20' * J;

stmC2S2 = stmC20 * stm0S2;

stmS1 = stt_i(:,:,1);
sttS1 = stt_i(:,:,:,1);

sttC20 = tensorCombine(stmS1, sttS1, stmC2S2, sttC2S2);
sttC10 = stt_t_i{1}(:,:,:,tcm_idx(1));
stt0C1 = invert_stt(stmC10, sttC10);

% finally combine them to get sttC2C1
sttC2C1 = tensorCombine(stm0C1, stt0C1, stmC20, sttC20);

Pc1_minus = P_i_minus(:,:,1);

R = simparams.R;
G = [zeros(3,3); eye(3,3)];
stmNC1 = stmN0 * stm0C1;
T1 = [-inv( stmNC1(1:3,4:6) ) * stmNC1(1:3,1:3), -eye(3)];

N1 = [zeros(3,6); T1];
IN1 = eye(6) + N1;


% Test dT1dx1_fd
dT1dx1_fd = zeros(3,6,6);

for i = 1:6
    x1dx = x(:);
    x1dx(i) = x(i) + dx;
    [stm_i1dx, stt_i1dx, x_i_f1dx, x_t1dx, stm_t1dx, stt_t_i1dx, t1dx, t_s1dx]  = createStateStmSttHistory(x1dx, simparams);
    [tcm_time1dx,tcm_idx1dx,min_tcm_dv1dx] = opt_multiple_tcm(x1dx, t1dx, stm_t1dx, simparams);
    [~, ~, tcm_dv_each1dx, P_i_minus1dx, P_i_plus1dx] = calc_covariance_tcmdv(x1dx, t1dx, stm_t1dx, tcm_time1dx, simparams);




    if simparams.target_final_maneuver
        final_seg = simparams.maneuverSegments(end);
        x1dx = reshape(x1dx,simparams.m,simparams.n);
        target_time = sum(x1dx(7,1:final_seg - 1));
        target_idx = find(target_time == t1dx) ;
        stmN01dx = stm_t1dx(:,:,target_idx);
        stm_t1dx = stm_t1dx(:,:,1:target_idx);
    else
        stmN01dx = stm_t1dx(:,:,end);
    end

    stmC10_1dx = stm_t1dx(:,:,tcm_idx1dx(1));
    stm0C1_1dx = -J * stmC10_1dx' * J;

    stmNC1_1dx = stmN01dx * stm0C1_1dx;

    T1_1dx = [-inv( stmNC1_1dx(1:3,4:6) ) * stmNC1_1dx(1:3,1:3), -eye(3)];

    dT1dx1_fd(:,:,i) = (T1_1dx - T1)./(dx);

    

end

%% having issues, testing sttN0
sttN0_fd = zeros(6,6,6);
for i = 1:6
    x1dx = x(:);
    x1dx(i) = x(i) + dx;
    [stm_i1dx, stt_i1dx, x_i_f1dx, x_t1dx, stm_t1dx, stt_t_i1dx, t1dx, t_s1dx]  = createStateStmSttHistory(x1dx, simparams);

    if simparams.target_final_maneuver
        final_seg = simparams.maneuverSegments(end);
        x1dx = reshape(x1dx,simparams.m,simparams.n);
        target_time = sum(x1dx(7,1:final_seg - 1));
        target_idx = find(target_time == t1dx) ;
        stmN01dx = stm_t1dx(:,:,target_idx);
        stm_t1dx = stm_t1dx(:,:,1:target_idx);
    else
        stmN01dx = stm_t1dx(:,:,end);
    end


    sttN0_fd = (stmN01dx - stmN0)./dx;

end
% stt20 = tensorCombine(stm_i(:,:,1),stt_i(:,:,:,1),stm_i(:,:,2),stt_i(:,:,:,2));
sttN0dx1 = tmult(stm_i(:,:,2),stt_i(:,:,:,1))


%%

% Get dT1dx1
% Need sttC10 and stt0C1, already have
% Need stmNC1f, the stm from the end of the correction segment to the
% target
seg1_t_idx = find((t_s==1)); % maybe a dangerous way to do this
stmC1f0 = stm_t(:,:,seg1_t_idx(end));
stm0C1f = -J * stmC1f0' * J;
stmNC1f = stmN0 * stm0C1f;

stti = stt_i(:,:,:,1);
dstmNC1dxi = tmult(stmNC1f, tmult(stti,stm0C1) + tmult(stm_i(:,:,1), stt0C1));
dT1dx1 = T_partial(stmNC1, dstmNC1dxi);
DIN1dx1 = [zeros(3,6,6); dT1dx1];




% calculate dstmC2C1dx1. the current situation is #2 [t0,i; tc1; tf,i; tc2]

stm20 = stm_t(:,:,find(t==x(7,1)));
% stm20 = stm_t(:,:,find(t==sum(x(7,1:2))));
stm02 = -J * stm20' * J;
stmC2ip1 = stmC20 * stm02;



dstmC2C1dx1 = tmult(  stmC2ip1, tmult(stm_i(:,:,1),stt0C1) + tmult(stt_i(:,:,:,1),stm0C1)  );

t1 = tmult(stmC2C1*G*R*G',dstmC2C1dx1,[0 1])                            + tmult(dstmC2C1dx1,G*R*G'*stmC2C1');

t2 = tmult(stmC2C1,tmult(IN1*Pc1_minus*IN1',dstmC2C1dx1,[0 1]))         + tmult(dstmC2C1dx1,tmult(IN1*Pc1_minus*IN1',stmC2C1'));

t3 = tmult(stmC2C1,tmult(DIN1dx1,Pc1_minus*IN1'*stmC2C1'))          + tmult(stmC2C1*IN1*Pc1_minus,tmult(DIN1dx1,stmC2C1,[1 1]));

t4 = tmult(stmC2C1*IN1,tmult(dPc1minusdx1_an,IN1'*stmC2C1'));


dPc2minusdx1_an = t1 + t2 + t3 + t4;

% in a function
dPc2minusdx1_anf = calc_dPckMinus(stmC2C1, dstmC2C1dx1, 2, simparams, T1, dT1dx1, Pc1_minus, dPc1minusdx1_an);




%% The above captured in the function sttCellCombine
% testing
% first test: if requested to match a segment, does it match
% [sttClast0, stmClast0] = sttCellCombine(stm_t, stt_t_i, t, t_s, 1, 90);
% stmClast0 - stm_i(:,:,1) % matches
% sttClast0 - stt_t_i{1}(:,:,:,end) % matches

% second test - second segment
seg_idx2 = find(t_s==2);
[sttClast0, stmClast0] = sttCellCombine(stm_t, stt_t_i, t, t_s, seg_idx2(1), seg_idx2(end));


stm_i2_0 = stm_t(:,:,seg_idx2(1)-1);
stm_0_i2 = -J * stm_i2_0' * J;
stm_i2f_0 = stm_t(:,:,seg_idx2(end));
stm_i2f_i2 = stm_i2f_0 * stm_0_i2;

stm_i2f_i2 - stm_i(:,:,2)   

stmClast0 - stm_i(:,:,2) % DISCREPANCY - WHAT IS GOING ON? DEBUG...%%%%% 10^-13 and 10^-14...but not zero
sttClast0 - stt_t_i{2}(:,:,:,end) % matches


% third test, segments 1 and 2
[sttClast0, stmClast0] = sttCellCombine(stm_t, stt_t_i, t, t_s, 1, seg_idx2(end));

stm21 = stm_i(:,:,2) * stm_i(:,:,1);
stmClast0-stm21
stmClast0 - stm_t(:,:,360)
stm21 - stm_t(:,:,360)

stt2 = stt_t_i{2}(:,:,:,end);
stt1 = stt_t_i{1}(:,:,:,end);
stm2 = stm_i(:,:,2);
stm1 = stm_i(:,:,1);

stt21 = tensorCombine(stm1,stt1,stm2,stt2);
sttClast0 - stt21


% fourth test, part of segment 3, segments 2 and 1
% [sttClast0, stmClast0] = sttCellCombine(stm_t, stt_t_i, t, t_s, 1, 400);
% 
% stt400_3 = stt_t_i{3}(:,:,:,41);
% stm400_0 = stm_t(:,:,400);
% stt0_2 = -J * stm21' * J;
% stm400_3 = stm400_0 * stt0_2;
% 
% stt400_0 = tensorCombine(stm21, stt21, stm400_3, stt400_3);
% 
% stt400_0 - sttClast0
% stm400_3*stm21 - stmClast0


%% testing the function calc_dstmCkClast
% [dstmCkClastdxi] = calc_dstmCkClast(x, t, t_s, stm_t, stm_i, stt_t_i, tCk, tClast, i)
% 

tClast = tcm_time(1);
% tCk = tcm_time(2);
tCk = t(80);

stmCk0 = stm_t(:,:,t==tCk);
stmClast0 = stm_t(:,:,find(t==tClast)-1);
stm0Clast = -J * stmClast0' * J;
stmCkClast = stmCk0 * stm0Clast;

dstmCkClast_fd = zeros(6,6,6);

for i = 7
    x1dx = x(:);
    x1dx(i) = x(i) + dx;
    [stm_i1dx, stt_i1dx, x_i_f1dx, x_t1dx, stm_t1dx, stt_t_i1dx, t1dx, t_s1dx]  = createStateStmSttHistory(x1dx, simparams);
    [x_t1dx, stm_t1dx, stt_t_i1dx, t1dx, t_s1dx] = addToStateStmSttHistory(x_t1dx, stm_t1dx, stt_t_i1dx, t1dx, t_s1dx, [tClast, tCk], simparams);
    %%%%%%%%%%%%%%%%%%%%%%%% AT SOME POINT, THE ABOVE FUNCTION NEEDS TO
    %%%%%%%%%%%%%%%%%%%%%%%% CORRECT TCM_IDX ALSO / INCORPORATE INTO THE
    %%%%%%%%%%%%%%%%%%%%%%%% UPDATE


    stmCk01dx = stm_t1dx(:,:,t1dx==tCk);
%     stmClast01dx = stm_t1dx(:,:,find(t1dx==tClast));
    stmClast01dx = stm_t1dx(:,:,find(t1dx==tClast)-1);
    stm0Clast1dx = -J * stmClast01dx' * J;
    stmCkClast1dx = stmCk01dx * stm0Clast1dx;

    dstmCkClast_fd(:,:,i) = (stmCkClast1dx - stmCkClast)./dx;

end



[dstmCkClast_an, dstmCkClastddti_an] = calc_dstmCkClast(x, x_i_f, t, t_s, stm_t, stm_i, stt_t_i, tCk, tClast, 1, simparams);
% Works for seg 1!!!!

%% Testing function calc_dstmNC


k = 2;

tCk = tcm_time(k);
tCk_idx = find(t==tCk);
tCk_seg = t_s(tCk_idx);



if simparams.target_final_maneuver
    final_seg = simparams.maneuverSegments(end);
    target_time = sum(x(7,1:final_seg - 1));
    target_idx = find(target_time == t) ;
    stmN0 = stm_t(:,:,target_idx);
    stm_t = stm_t(:,:,1:target_idx);
    t = t(1:target_idx);
    t_s = t_s(1:target_idx);
else
    stmN0 = stm_t(:,:,end);
end


stmCk0 = stm_t(:,:,tCk_idx);
stm0Ck = -J * stmCk0' * J;

stmNCk = stmN0 * stm0Ck;

for i = 1:21
    x1dx = x(:);
    x1dx(i) = x(i) + dx;
    [stm_i1dx, stt_i1dx, x_i_f1dx, x_t1dx, stm_t1dx, stt_t_i1dx, t1dx, t_s1dx]  = createStateStmSttHistory(x1dx, simparams);
    [x_t1dx, stm_t1dx, stt_t_i1dx, t1dx, t_s1dx] = addToStateStmSttHistory(x_t1dx, stm_t1dx, stt_t_i1dx, t1dx, t_s1dx, tcm_time, simparams);


    x1dx  = reshape(x1dx,simparams.m,simparams.n);

    if simparams.target_final_maneuver
        final_seg = simparams.maneuverSegments(end);
        target_time1dx = sum(x1dx(7,1:final_seg - 1));
        target_idx1dx = find(target_time1dx == t1dx);
        stmN01dx = stm_t1dx(:,:,target_idx1dx);
        stm_t1dx = stm_t1dx(:,:,1:target_idx1dx);
        t1dx = t1dx(1:target_idx1dx);
        t_s1dx = t_s1dx(1:target_idx1dx);
    else
        stmN01dx = stm_t1dx(:,:,end);
    end




    tCk1dx = tcm_time(k);
    tCk_idx1dx = find(t1dx==tCk1dx);
    x1dx = reshape(x1dx,simparams.m, simparams.n);
%     target_time1dx = sum(x1dx(7,1:final_seg - 1));


    stmCk01dx = stm_t1dx(:,:,tCk_idx1dx);
    stm0Ck1dx = -J * stmCk01dx' * J;

    stmNCk1dx = stmN01dx * stm0Ck1dx;

    dstmNC_fd(:,:,i) = (stmNCk1dx - stmNCk)./dx;

    %%%%%%%%%%% second approach to finite difference test
    % don't change state or STM propagation from t0 to 1
    % propagate from t1 to tCk
    % propagate from tCk to tN

end
[dstmNCdxi,dstmNCddti] = calc_dstmNC(x, x_i_f, t, t_s, stm_t, stm_i, stt_t_i, tCk, tCk_idx, tCk_seg, stmN0, 2, simparams);


[tcm_gradient] = calc_multiple_tcm_gradient(x, x_t, x_i_f, stm_i, stt_i, stm_t, stm_t_i, stt_t_i, t, t_s, tcm_time, tcm_idx, tcm_dv_each, P_i_minus, P_i_plus, simparams);

%% second appraoch at verifying calc_dstmNC
    %%%%%%%%%%% second approach to finite difference test
    % don't change state or STM propagation from t0 to 1
    % propagate from t1 to tCk
    % propagate from tCk to tN

%


k = 2;

tCk = tcm_time(k);
tCk_idx = find(t==tCk);
tCk_seg = t_s(tCk_idx);

simpar2 = simparams;
simpar2.n=1;
[stm_i1dx, stt_i1dx, x_i_f1dx, x_t1dx, stm_t1dx, stt_t_i1dx, t1dx, t_s1dx]  = createStateStmSttHistory(x(1:7), simpar2);
stm1f0 = stm_t1dx(:,:,end);

x20 = x(1:6,2);

% starting from the end, propagate dT from t1,0 to tCk

dt = tCk - x(7);

[stm_i1dx, stt_i1dx, x_i_f1dx, x_t1dx, stm_t1dx, stt_t_i1dx, t1dx, t_s1dx]  = createStateStmSttHistory([x20; dt], simpar2);

stmCk1f = stm_t1dx(:,:,end);
xCk = x_i_f1dx;
dt2 = x(14) - dt;

% propagate the rest of the way to 2f


[stm_i1dx, stt_i1dx, x_i_f1dx, x_t1dx, stm_t1dx, stt_t_i1dx, t1dx, t_s1dx]  = createStateStmSttHistory([xCk; dt2], simpar2);

stmNCk_t = stm_i1dx;





for i = 1:6
    x20 = x(1:6,2);
    x20(i) = x20(i) + dx;
    
    % starting from the end, propagate dT from t1,0 to tCk
    
    dt = tCk - x(7);
    
    [stm_i1dx, stt_i1dx, x_i_f1dx, x_t1dx, stm_t1dx, stt_t_i1dx, t1dx, t_s1dx]  = createStateStmSttHistory([x20; dt], simpar2);
    
    stmCk1f = stm_t1dx(:,:,end);
    xCk = x_i_f1dx;
    dt2 = x(14) - dt;
    
    % propagate the rest of the way to 2f
    
    
    [stm_i1dx, stt_i1dx, x_i_f1dx, x_t1dx, stm_t1dx, stt_t_i1dx, t1dx, t_s1dx]  = createStateStmSttHistory([xCk; dt2], simpar2);
    
    stmNCk_tdx = stm_i1dx;
    dstmNCktdx_fd(:,:,i) = (stmNCk_tdx - stmNCk_t)./dx;
end

% Gets the same result - it appears the problem isn't with the test...

% assert(0)

%% try to verify the tcm r gradient

stmC10 = stm_t(:,:,tcm_idx(1));
sttC10 = stt_t_i{1}(:,:,:,tcm_idx(1));



seg2_c2_idx = tcm_idx(2) - size(stt_t_i{1},4) + 1; %%%%% THIS COULD BE ONE OFF
sttC2S2 = stt_t_i{t_s(tcm_idx(2))}(:,:,:,seg2_c2_idx);

seg2_indices = find(t_s==2);
seg2_start_idx = seg2_indices(1);
stmS20 = stm_t(:,:,seg2_start_idx);
stm0S2 = -J * stmS20' * J;

stmC2S2 = stmC20 * stm0S2;

stmS1 = stt_i(:,:,1);
sttS1 = stt_i(:,:,:,1);




sttC20 = tensorCombine(stmS1, sttS1, stmC2S2, sttC2S2);
sttC10 = stt_t_i{1}(:,:,:,tcm_idx(1));
stt0C1 = invert_stt(stmC10, sttC10);

% finally combine them to get sttC2C1
sttC2C1 = tensorCombine(stm0C1, stt0C1, stmC20, sttC20);

Pc1_minus = P_i_minus(:,:,1);

R = simparams.R;
G = [zeros(3,3); eye(3,3)];
stmNC1 = stmN0 * stm0C1;
T1 = [-inv( stmNC1(1:3,4:6) ) * stmNC1(1:3,1:3), -eye(3)];

N1 = [zeros(3,6); T1];
IN1 = eye(6) + N1;

if simparams.target_final_maneuver
    final_seg = simparams.maneuverSegments(end);
    target_time = sum(x(7,1:final_seg - 1));
    target_idx = find(target_time == t) ;
    stmN0 = stm_t(:,:,target_idx);
    stm_t = stm_t(:,:,1:target_idx);
else
    stmN0 = stm_t(:,:,end);
end




k=2;
% stmCk0 = stm_t(:,:,tcm_idx(k));
% stm0Ck = -J * stmCk0' * J;
% stmNCk = stmN0 * stm0Ck;
target_idx=360;
stmNCk = dynCellCombine(t, t_s, tcm_idx(k), target_idx, stm_t_i);


Tk = [-inv( stmNCk(1:3,4:6) ) * stmNCk(1:3,1:3), -eye(3)];
Pck_minus = P_i_minus(:,:,k);
sig_dvk = sqrt(trace(Tk * Pck_minus * Tk'));

[~, tcm_dv_total, tcm_dv] = calc_covariance_tcmdv(x, t, stm_t, tcm_time, simparams);
tcmR_dv_total = sum(tcm_dv(1:end-1));
tcm1_dv = tcm_dv(1);
tcm2_dv = tcm_dv(2);
tcm3_dv = tcm_dv(3);
dtcmR_fd = zeros(length(x(:)),1);
dtcmV_fd = zeros(length(x(:)),1);
dtcm_total_fd=zeros(length(x(:)),1);
dTkdxi = zeros(3,6,length(x(:)));
dstmNCk_fd = zeros(6,6,length(x(:)));


t20_idx = find(t==x(7,1));
t2f_idx = find(t==sum(x(7,1:2)));
[stm_20_tc2, stt_20_tc2]  = dynCellCombine(t, t_s, tcm_idx(k), t20_idx, stm_t_i, stt_t_i);
stm_2f_20  = dynCellCombine(t, t_s, t20_idx, t2f_idx, stm_t_i);



if k == 1
    tClast_idx = 1;
else
    tClast_idx = find(t==tcm_time(k-1));
end
tCk_idx = find(t==tcm_time(k));
stmCkClast = dynCellCombine(t, t_s, tClast_idx, tCk_idx, stm_t_i);



for i = 1:21
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
    dtcmR_fd(i) = (tcmR_dv_total1dx-tcmR_dv_total)/dx;
    dtcmV_fd(i) = (tcm_dv1dx(end)-tcm_dv(end))/dx;

    dtcm_total_fd(i) = (tcm_dv_total1dx - tcm_dv_total)./dx;

    



%     stmCk0_1dx = stm_t1dx(:,:,t1dx==tcm_time(k));
%     stm0Ck_1dx = - J * stmCk0_1dx' * J;
% 
%     stmNCk1dx = stmN01dx * stm0Ck_1dx; % this is including mods to stmN0 from segment 1, when the correction is in seg 2, for example, which is incorrect
    tCk_idx = find(t1dx==tcm_time(k));
    stmNCk1dx = dynCellCombine(t1dx, t_s1dx, tCk_idx, target_idx1dx, stm_t_i1dx);

    if k == 1
        tClast_idx = 1;
    else
        tClast_idx = find(t==tcm_time(k-1));
    end

    

    dstmNCk_fd(:,:,i) = (stmNCk1dx - stmNCk)./dx;

    Tk1dx = [-inv( stmNCk1dx(1:3,4:6) ) * stmNCk1dx(1:3,1:3), -eye(3)];

    Pck_minus1dx = P_i_minus1dx(:,:,k);
    sig_dvk1dx = sqrt(trace(Tk1dx * Pck_minus1dx * Tk1dx'));

    dsig_dvk_dx1dx_fd(i) = (sig_dvk1dx - sig_dvk)./dx;

    dTkdxi(:,:,i) = (Tk1dx - Tk)./dx;




    t20_idx1dx = find(t1dx==x1dx(7,1));
    stm_20_tc21dx  = dynCellCombine(t1dx, t_s1dx, tcm_idx1dx(k), t20_idx1dx, stm_t_i1dx);

    t2f_idx1dx = find(t1dx==sum(x1dx(7,1:2)));
    stm_2f_201dx  = dynCellCombine(t1dx, t_s1dx, t20_idx1dx, t2f_idx1dx, stm_t_i1dx);



    dstm_20_tc2_fd(:,:,i) = (stm_20_tc21dx - stm_20_tc2)./dx;
    dstm_2f_20_fd(:,:,i) = (stm_2f_201dx - stm_2f_20)./dx;



    % testing dstmCkClast_dxi


    if k == 1
        tClast_idx1dx = 1;
    else
        tClast_idx1dx = find(t1dx==tcm_time1dx(k-1));
    end
    tCk_idx1dx = find(t1dx==tcm_time1dx(k));
    stmCkClast1dx = dynCellCombine(t1dx, t_s1dx, tClast_idx1dx, tCk_idx1dx, stm_t_i1dx);


    dstmCkClast_dxi(:,:,i) = (stmCkClast1dx - stmCkClast)./dx;


end

dstm_20_tc2_dx2 = dstm_20_tc2_fd(:,:,8:13); %%%%% bug in dynCellCombine getting the same result MAYBE!!!!!!
dstm_2f_20_dx2 = dstm_2f_20_fd(:,:,8:13);
dstmNC2dx2_verify_fd = tmult(dstm_2f_20_dx2, stm_20_tc2) + tmult(stm_2f_20, dstm_20_tc2_dx2);



% trying to get the same dstm_20_tc2_dx2 as above
tc2_cell_idx = find_cell_t_i_idx(t,t_s,stt_t_i,tcm_idx(2));
stt_c2_t0 = stt_t_i{2}(:,:,:,tc2_cell_idx);
stm_c2_t0 = stm_t_i{2}(:,:,tc2_cell_idx);
stm_t0_c2 = -J * stm_c2_t0' * J;

stt_t0_c2 = invert_stt(stm_c2_t0,stt_c2_t0)

% test inverse matrix property
dstmc2_20_inv_dx2 = -tmult(stm_t0_c2,  tmult(stt_c2_t0, stm_t0_c2)  )
% it worked...but why did the other way work for dx1 when k=1????
%%%% BIG TODO: CLEAN UP IMPLEMENTATIONS...AND TEST...AND FIGURE OUT WHY THE
%%%% OTHER WAY WORKED THE FIRST TIME


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

[dtcm1_fd' ,  dtcm2_fd' ,  dtcm3_fd']


% easier to debug into the calc_multiple_tcm_gradient function and directly
% observe the tcm_gradient_r values
[tcm_gradient] = calc_multiple_tcm_gradient(x, x_t, x_i_f, stm_i, stt_i, stm_t, stm_t_i, stt_t_i, t, t_s, tcm_time, tcm_idx, tcm_dv_each, P_i_minus, P_i_plus, simparams);
