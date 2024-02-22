function [xfinal, stm, stt, xstmstt_t, ti] = stateStmSttProp(x_initial, t, simparams,  stm_initial, stt_initial)

mu = simparams.mu;
options = simparams.options;
dynSys = simparams.dynSys;
nsv = simparams.nsv; % Number of state variables


if ~( strcmp(dynSys,'2bp') || strcmp(dynSys,'cr3bp') || strcmp(dynSys,'br4bp_em') || strcmp(dynSys,'br4bp_sb1'))
    error('Error: Dynamical system identifier text not recognized!')
end


if nargin == 3
    stm_initial = eye(nsv);
    stt_initial = zeros(nsv, nsv, nsv);
end

if nargin == 4
    stt_initial = zeros(nsv, nsv, nsv);
end

xstmstt_initial = [x_initial; stm_initial(:); stt_initial(:)];
if t == 0
    xfinal = x_initial;
    stm = stm_initial;
    stt = stt_initial;
    xstmstt_t = xstmstt_initial;
else
    if strcmp(dynSys,'2bp')
    % Two body numerical propagation of state, STM, and STT
        [ti, xstmstt_t] = ode113(@r2bp_stt2_de, [0,t], xstmstt_initial, options, mu);

    elseif strcmp(dynSys,'cr3bp')
        [ti, xstmstt_t] = ode113(@cr3bp_sFrame_nd_stt2_de, [0,t], xstmstt_initial, options, mu);

    elseif strcmp(dynSys,'br4bp_em')
        % Bicircular restricted four body problem, earth-moon frame
        assert(0,'STT propagation in BR4BP EM frame not implemented yet.')

    elseif strcmp(dynSys,'br4bp_sb1')
        % Bicircular restricted four body problem, Sun-B1 frame        
        [ti, xstmstt_t] = ode113(@br4bp_sb1_state_stm_stt_de, [0,t], xstmstt_initial, options, simparams.mub, mu, simparams.a4, simparams.theta_em_dot);

    end

    xstmstt_f = xstmstt_t(end,:);

    xfinal = reshape(xstmstt_f(1:nsv),nsv,1);
    stm = reshape(xstmstt_f(nsv+1:nsv+nsv^2),nsv,nsv);
    stt = reshape(xstmstt_f(nsv+nsv^2+1:nsv+nsv^2+nsv^3),nsv,nsv,nsv);
end






end