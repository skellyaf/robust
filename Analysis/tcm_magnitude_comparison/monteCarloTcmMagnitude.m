function [tcm_dv_mag] = monteCarloTcmMagnitude(x0_i, x_target, t, tcm_time, stm_t, P0, num_samples, percentile, simparams)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

options = optimset('Display','off');
options.MaxFunEvals = 3e3;
% options.MaxIterations = 1e3;

total_time = t(end);
tcm_idx = find(t==tcm_time);


stmN0 = stm_t(:,:,end);
stmC0 = stm_t(:,:,tcm_idx);
stm0C = invert_stm(stmC0, simparams);
stmNC = stmN0 * stm0C;


T = [-inv( stmNC(1:3,4:6) ) * stmNC(1:3,1:3), -eye(3)];

tcm_mag_samples = zeros(1,num_samples);

x_tcm = stateProp(x0_i, tcm_time, simparams);

x_final = stateProp(x_tcm, total_time - tcm_time, simparams);


%% covariance matrix, rotation, and eigenvectors
P0_i = P0;
[eigvec, P0_p] = eig(P0_i);

T_p2i = eigvec;
T_i2p = eigvec';


%% samples


for i = 1:num_samples
    x0_p = T_i2p * x0_i;

    x0r_p = x0_p + randn(6,1).*sqrt(diag(P0_p));
    x0r_i = T_p2i * x0r_p;

%     x0r = x0 + randn(6,1).*sqrt(diag(P0));
    xr_tcm = stateProp(x0r_i, tcm_time, simparams);
    delx = xr_tcm - x_tcm;
    tcm_dv = T * delx;

    xr_c_guess = xr_tcm + [zeros(3,1); tcm_dv];

%     xr_notcm_final = stateProp(xr_tcm, total_time - tcm_time, simparams);


    xr_final = stateProp(xr_c_guess, total_time - tcm_time, simparams);

%     dx = xr_final - x_final
%     dxnotcm = xr_notcm_final - x_final


    [tcm_dv_fsolve, fval, exitflag, output] = fsolve( @(tcm_dv)fsolve_tcm_target_function2(tcm_dv, xr_tcm, x_final, total_time - tcm_time, simparams), tcm_dv, options);

    if exitflag < 1
        ppp=1;
        tcm_mag_samples(i) = tcm_mag_samples(i-1);
    end

    tcm_mag_samples(i) = vecnorm(tcm_dv_fsolve);



end

tcm_dv_mag = prctile(tcm_mag_samples, percentile);