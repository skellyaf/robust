function [xfinal, stm, xstm_t, ti] = stateStmProp(xinitial, t, simparams,  stm_initial)

mu = simparams.mu;
options = simparams.options;
dynSys = simparams.dynSys;

if ~( strcmp(dynSys,'2bp') || strcmp(dynSys,'cr3bp') )
    error('Error: Dynamical system identifier text not recognized!')
end


if nargin == 3
    stm_initial = eye(6);
end
if length(stm_initial(:))~=36
    assert(0);
end
xstm = [xinitial; reshape(stm_initial, 36, 1)];
if t ~= 0
    if strcmp(dynSys,'2bp')
    % Two body numerical propagation of state and STM
        [ti, xstm_t] = ode113(@r2bp_stm_de, [0,t], xstm, options, mu);
    elseif strcmp(dynSys,'cr3bp')
    % Circular restricted three body numerical propagation of state and STM
        [ti, xstm_t] = ode113(@cr3bp_sFrame_nd_stm_de, [0,t], xstm, options, mu);
    end
    stm = reshape(xstm_t(end,7:42),6,6);    
    xfinal = reshape(xstm_t(end,1:6),6,1);
else
    stm = stm_initial;
    xfinal = xinitial;
    xstm_t = xstm;
end
end