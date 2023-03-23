function [xfinal, x_t, t] = stateProp(xinitial, t, simparams)

options = simparams.options;
dynSys = simparams.dynSys;
mu = simparams.mu;

if ~( strcmp(dynSys,'2bp') || strcmp(dynSys,'cr3bp') )
    error('Error: Dynamical system identifier text not recognized!')
end


if t ~= 0

    if strcmp(dynSys,'2bp')
        % Two body numerical propagation of state
        [t, x_t] = ode113(@r2bp_de, [0,t], xinitial, options, mu);
    elseif strcmp(dynSys,'cr3bp')
        % Circular restricted three body numerical propagation of state
        [t, x_t] = ode113(@cr3bp_sFrame_nd_de, [0,t], xinitial, options, mu);
    end

    xfinal = reshape(x_t(end,1:6),6,1);
else
    xfinal = xinitial;
    x_t = xinitial;
end
end