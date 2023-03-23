function [tcm_total, tcm_time] = concurrentManeuverTcm(x, t, stm_t, simparams)
%minTcmPair_rv Calculates the TCM 3 sigma magnitude at a specific point
%along a trajectory, concurrent with one of the maneuvers

P_initial = simparams.P_initial;
x = reshape(x,simparams.m,simparams.n);





stmN0 = stm_t(:,:,end);

% Symplectic unit matrix
J = [zeros(3,3), eye(3,3); -eye(3,3), zeros(3,3)];




%     t_c = t(i);
%     x_c = x_t(i,:)';
stmC0 = stm_t(:,:,tcm_idx);
stm0C = -J * stmC0' * J; % symplectic inverse
stmNC = stmN0 * stm0C;

T = [-inv( stmNC(1:3,4:6) ) * stmNC(1:3,1:3), -eye(3)];

% Magnitude of correction to final position dispersion
dvR3sigma = 3 * sqrt( trace( T * stmC0 * P_initial * stmC0' * T' ) );

% Magnitude of correction to final velocity dispersion
W = [stmNC(4:6,4:6), -stmNC(4:6,1:3)];
L = [W*T', zeros(3,3)];
dvV3sigma = 3 * sqrt( trace( L * stmC0 * P_initial * stmC0' * L' ) );


tcm_total = dvR3sigma + dvV3sigma;


end