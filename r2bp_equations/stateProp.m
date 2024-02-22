function [xfinal, x_t, t] = stateProp(xinitial, t, simparams)

options = simparams.options;
dynSys = simparams.dynSys;
mu = simparams.mu;
nsv = simparams.m - 1; % number of state variables

if ~( strcmp(dynSys,'2bp') || strcmp(dynSys,'cr3bp') || strcmp(dynSys,'br4bp_sb1'))
    error('Error: Dynamical system identifier text not recognized!')
end


if t ~= 0

    if strcmp(dynSys,'2bp')
        % Two body numerical propagation of state
        [t, x_t] = ode113(@r2bp_de, [0,t], xinitial, options, mu);
    elseif strcmp(dynSys,'cr3bp')
        % Circular restricted three body numerical propagation of state
        [t, x_t] = ode113(@cr3bp_sFrame_nd_de, [0,t], xinitial, options, mu);
    elseif strcmp(dynSys,'br4bp_sb1')
        % Bicircular restricted four body problem, Sun-B1 frame
        [t, x_t] = ode113(@br4bp_sb1_state_de, [0,t], xinitial, options, simparams.mub, mu, simparams.a4, simparams.theta_em_dot);

    end

    xfinal = reshape(x_t(end,1:nsv),nsv,1);
else
    xfinal = xinitial;
    x_t = xinitial;
end
end