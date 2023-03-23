function [tcm_min, tcm_time] = minTcmPair_rv(t, stm_t, simparams)
%minTcmPair_rv Calculates the minimum TCM sum along a trajectory
% There are two TCMs required to hit a target state. The first occurs along
% the trajectory at an optimal execution time which targets a zero final
% position dispersion. The second TCM occurs at the end of the trajectory
% and corrects any final velocity dispersion.

P_initial = simparams.P_initial;
x = reshape(x,simparams.m,simparams.n);




% Propagate entire trajectory once to get stmNC
stmN0 = stm_t(:,:,end);

% Symplectic unit matrix
J = [zeros(3,3), eye(3,3); -eye(3,3), zeros(3,3)];

% Initialization
dvR3sigma = zeros(size(t,1),1);
dvV3sigma = zeros(size(t,1),1);

%%%%% there may be a bug here-----no dv1dvctied param exists anymore!!!!!
if simparams.dv1dvctied
    % Only run through the for loop below once, at the correct index for
    % when the first maneuver occurs
    idx1 = %%%%%%%%%%%%%%%%% WHAT HAPPENED HERE???????????
    idx2 = idx1;
else
    % Otherwise, run through the entire trajectory to find the lowest TCM
    idx1 = 1;
    idx2 = size(t,1);
end

% Calculate the minimum TCM pair
for i = 1:size(t,1)

%     t_c = t(i);
%     x_c = x_t(i,:)';
    stmC0 = stm_t(:,:,i);
    stm0C = -J * stmC0' * J; % symplectic inverse
    stmNC = stmN0 * stm0C;

    T = [-inv( stmNC(1:3,4:6) ) * stmNC(1:3,1:3), -eye(3)];

    % Magnitude of correction to final position dispersion
    dvR3sigma(i) = 3 * sqrt( trace( T * stmC0 * P_initial * stmC0' * T' ) );

    % Magnitude of correction to final velocity dispersion
    W = [stmNC(4:6,4:6), -stmNC(4:6,1:3)];
    L = [W*T', zeros(3,3)];
    dvV3sigma(i) = 3 * sqrt( trace( L * stmC0 * P_initial * stmC0' * L' ) );

end

tcm_total = dvR3sigma + dvV3sigma;
[tcm_min, tcm_min_idx] = min(tcm_total);

tcm_time = t(tcm_min_idx);



end