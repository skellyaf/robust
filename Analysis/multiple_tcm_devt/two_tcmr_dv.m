function [two_tcm_total_dv,tcm1_dv, tcm2_dv, tcm3_dv] = two_tcmr_dv(x, t, stm_t, tcm_time, simparams)
%two_tcmr_dv computes the total delta V RSS along a reference trajectory
%given the times the two TCMs occur (plus the final velocity correction
%delta V. The returned values are the RSS of the TCM covariance matrix,
%which is an estimate for the standard deviation of the magnitude of the
%TCM.

assert(length(tcm_time) == 2,'Wrong number of TCM times (not 2) pass to function');

% Symplectic unit matrix
J = [zeros(3,3), eye(3,3); -eye(3,3), zeros(3,3)];
P_i = simparams.P_initial;


% calc stmN0
if simparams.target_final_maneuver
    final_seg = simparams.maneuverSegments(end);
    target_time = sum(x(7,1:final_seg - 1));
    target_idx = find(target_time == t) ;
    stmN0 = stm_t(:,:,target_idx);

else
    stmN0 = stm_t(:,:,end);
end


[P_target, P_t] = calc_covariance_history(x, t, stm_t, tcm_time, simparams); % Can make faster by not calculating entire history...likely unnecessary here
t1_idx = find(t==tcm_time(1));
stmC10 = stm_t(:,:,t1_idx); 
stm0C1 = -J * stmC10' * J; % symplectic inverse 
stmNC1 = stmN0 * stm0C1;

T1 = [-inv( stmNC1(1:3,4:6) ) * stmNC1(1:3,1:3), -eye(3)];
P_tcm1 = T1 * stmC10 * P_i * stmC10' * T1' + simparams.R;
tcm1_dv = sqrt(trace(P_tcm1));




t2_idx = find(t==tcm_time(2));
stmC20 = stm_t(:,:,t2_idx); 
stm0C2 = -J * stmC20' * J; % symplectic inverse 
stmNC2 = stmN0 * stm0C2;

T2 = [-inv( stmNC2(1:3,4:6) ) * stmNC2(1:3,1:3), -eye(3)];

stmC2C1 = stmC20 * stm0C1;

P_tcm2 = T2 * stmC2C1 * P_t(:,:,t1_idx) * stmC2C1' * T2' + simparams.R; 
tcm2_dv = sqrt(trace(P_tcm2));

tcm3_dv = sqrt(trace(P_t(4:6,4:6,end)));

two_tcm_total_dv = tcm1_dv + tcm2_dv + tcm3_dv;

end