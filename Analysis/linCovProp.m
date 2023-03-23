function [P_t] = linCovProp(x, stm_t, t, tcm_time, simparams)
%linCovProp uses linear covariance propagation to propagate an initial
%state dispersion throughout a trajectory defined by stm_t

%%%% TO DO: THE stm_i AND stm_t ARRAYS OF STM'S ARE SET UP DIFFERENTLY.
%%%% NEED TO DO SOMETHING TO ACCOMMODATE HERE, EITHER ONE, OR CHANGE HOW
%%%% THEY ARE MADE. EITHER WAY - SHOULD PROBABLY BE STANDARDIZED ACROSS THE
%%%% CODE.

% Reshape x
x = reshape(x,simparams.m,simparams.n);

% The index of t and stm_t where the TCM is performed
tcm_t_idx = find(tcm_time == t);

% Symplectic unit matrix
J = [zeros(3,3), eye(3,3); -eye(3,3), zeros(3,3)];

stmC0 = stm_t(:,:,tcm_t_idx);
stm0C = -J * stmC0' * J; % symplectic inverse

if simparams.target_final_maneuver
    % Index of final maneuver
    final_dv_idx = simparams.maneuverSegments(end); % the seg w/an impulsive maneuver at its beginning
    % Find the time that the final maneuver occurs
    t_final_dv = sum(x(7,1:final_dv_idx-1));
    % Find the corresponding t index of the final maneuver time
    final_dv_t_idx = find(t_final_dv == t);
    % The "end"/target of the correction is then based on the target index
    stmN0 = stm_t(:,:,final_dv_t_idx);
else
    
    % If the target is the end of the trajectory, the end of the stm_t
    % tensor is used:
    stmN0 = stm_t(:,:,end);
end

% The dynamics from the correction to the target are then:
stmNC = stmN0 * stm0C;

% The STM targeting matrix is
if simparams.perform_correction
    T = [-inv( stmNC(1:3,4:6) ) * stmNC(1:3,1:3), -eye(3)];
    N = [zeros(3,6); T];
    IN = eye(6) + N;
else
    IN = eye(6);
end

%% Vectorizing the propagation rather than using a for loop

stm_t_c0 = stm_t(:,:,1:tcm_t_idx);
P_c0 =  tmult(stm_t_c0, tmult(simparams.P_initial, stm_t_c0, [0 1]));
P_c = tmult(IN, tmult(P_c0(:,:,end), IN, [0 1]));

stm_t_nc = tmult(stm_t(:,:,tcm_t_idx+1:final_dv_t_idx), stm0C);
P_nc = tmult(stm_t_nc, tmult(P_c, stm_t_nc, [0 1]));

P_n = P_nc(:,:,end);
P_n(4:6,4:6) = zeros(3,3);

if size(stm_t,3) > final_dv_t_idx

    stm0N = -J * stmN0' * J;
    stm_t_fn = tmult(stm_t(:,:,final_dv_t_idx+1:end), stm0N);
    P_fn = tmult(stm_t_fn, tmult(P_n, stm_t_fn, [0 1]));
    
    P_t = cat(3,P_c0, P_nc, P_fn);
else
    P_t = cat(3,P_c0, P_nc);
end





end