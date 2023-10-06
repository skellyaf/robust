function [xfinal, stm, Q, xstmQ_t, ti] = statestmQProp(xinitial, t, simparams, stm_initial, Q_initial)

mu = simparams.mu;
Qt = simparams.Qt;
options = simparams.options;
dynSys = simparams.dynSys;

if ~( strcmp(dynSys,'2bp') || strcmp(dynSys,'cr3bp') )
    error('Error: Dynamical system identifier text not recognized!')
end


if nargin == 3
    stm_initial = eye(6);
    Q_initial = zeros(6,6);
end
if nargin == 4
    Q_initial = zeros(6,6);
end


if length(Q_initial(:))~=36
    assert(0);
end
xstmQ = [xinitial; stm_initial(:); Q_initial(:)];
if t ~= 0
    if strcmp(dynSys,'2bp')
    % Two body numerical propagation of state and STM
        assert(0,'not set up yet!');
        [ti, xstmQ_t] = ode113(@r2bp_stm_de, [0,t], xQQ2, options, mu);
    elseif strcmp(dynSys,'cr3bp')
    % Circular restricted three body numerical propagation of state and STM
        [ti, xstmQ_t] = ode113(@statestmQdot, [0,t], xstmQ, options, mu, Qt);
    end

    xfinal = reshape(xstmQ_t(end,1:6),6,1);
    stm = reshape(xstmQ_t(end,7:42),6,6);
    Q = reshape(xstmQ_t(end,43:78),6,6);
    
else
    Q = Q_initial;
    xfinal = xinitial;
    xstmQ_t = xstmQ;
end
end