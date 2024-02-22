function [xfinal, stm, Q, xstmQ_t, ti] = statestmQProp(xinitial, t, simparams, stm_initial, Q_initial)

mu = simparams.mu;
Qt = simparams.Qt;
options = simparams.options;
dynSys = simparams.dynSys;
nsv = simparams.nsv;

if ~( strcmp(dynSys,'2bp') || strcmp(dynSys,'cr3bp') || strcmp(dynSys,'br4bp_em') || strcmp(dynSys,'br4bp_sb1') )
    error('Error: Dynamical system identifier text not recognized!')
end


if nargin == 3
    stm_initial = eye(nsv);
    Q_initial = zeros(nsv,nsv);
end
if nargin == 4
    Q_initial = zeros(nsv,nsv);
end


if length(Q_initial(:))~=nsv^2
    assert(0);
end
xstmQ = [xinitial; stm_initial(:); Q_initial(:)];

if t ~= 0
    if strcmp(dynSys,'2bp')
    % Two body numerical propagation of state and STM
        assert(0,'not set up yet!');
%         [ti, xstmQ_t] = ode113(@r2bp_stm_de, [0,t], xQQ2, options, mu);

    elseif strcmp(dynSys,'cr3bp')
    % Circular restricted three body numerical propagation of state and STM
        [ti, xstmQ_t] = ode113(@statestmQdot, [0,t], xstmQ, options, mu, Qt);

    elseif strcmp(dynSys,'br4bp_em')
        % Bicircular restricted four body problem, earth-moon frame
        assert(0,'not set up yet!');
%         [ti, xstm_t] = ode113(@br4bp_em_state_stm_de, [0,t], xstm, options, mu, simparams.m4, simparams.a4, simparams.theta_s_dot);

    elseif strcmp(dynSys,'br4bp_sb1')
        % Bicircular restricted four body problem, Sun-B1 frame
        [ti, xstmQ_t] = ode113(@br4bp_sb1_state_stm_Q_de, [0,t], xstmQ, options, simparams.mub, mu, simparams.a4, simparams.theta_em_dot, Qt);

    end

    xfinal = reshape(xstmQ_t(end,1:nsv),nsv,1);
    stm = reshape(xstmQ_t(end,nsv+1:nsv+nsv*nsv),nsv,nsv);
    Q = reshape(xstmQ_t(end,nsv+nsv*nsv+1:nsv+nsv*nsv+nsv*nsv),nsv,nsv);
    
else
    Q = Q_initial;
    xfinal = xinitial;
    xstmQ_t = xstmQ;
end
end