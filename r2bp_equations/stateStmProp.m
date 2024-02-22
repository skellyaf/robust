function [xfinal, stm, xstm_t, ti] = stateStmProp(xinitial, t, simparams,  stm_initial)

mu = simparams.mu;
options = simparams.options;
dynSys = simparams.dynSys;
nsv = simparams.m - 1; % number of state variables

if ~( strcmp(dynSys,'2bp') || strcmp(dynSys,'cr3bp') || strcmp(dynSys,'br4bp_em') || strcmp(dynSys,'br4bp_sb1'))
    error('Error: Dynamical system identifier text not recognized!')
end


if nargin == 3
    stm_initial = eye(nsv);
end

if length(stm_initial(:)) ~= (nsv)^2
    assert(0);
end


xstm = [xinitial; stm_initial(:)];


if t ~= 0
    if strcmp(dynSys,'2bp')
        % Two body numerical propagation of state and STM
        [ti, xstm_t] = ode113(@r2bp_stm_de, [0,t], xstm, options, mu);

    elseif strcmp(dynSys,'cr3bp')
        % Circular restricted three body numerical propagation of state and STM
        [ti, xstm_t] = ode113(@cr3bp_sFrame_nd_stm_de, [0,t], xstm, options, mu);

    elseif strcmp(dynSys,'br4bp_em')
        % Bicircular restricted four body problem, earth-moon frame
        [ti, xstm_t] = ode113(@br4bp_em_state_stm_de, [0,t], xstm, options, mu, simparams.m4, simparams.a4, simparams.theta_s_dot);

    elseif strcmp(dynSys,'br4bp_sb1')
        % Bicircular restricted four body problem, Sun-B1 frame
        [ti, xstm_t] = ode113(@br4bp_sb1_state_stm_de, [0,t], xstm, options, simparams.mub, mu, simparams.a4, simparams.theta_em_dot);

    end

    xfinal = reshape(xstm_t(end,1:nsv),nsv,1);
    stm = reshape(xstm_t(end,nsv+1:nsv+nsv^2),nsv,nsv);
    
else
    stm = stm_initial;
    xfinal = xinitial;
    xstm_t = xstm';
end
end