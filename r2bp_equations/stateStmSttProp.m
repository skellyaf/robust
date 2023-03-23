function [xfinal, stm, stt, xstmstt_t, ti] = stateStmSttProp(x_initial, t, simparams,  stm_initial, stt_initial)

mu = simparams.mu;
options = simparams.options;
dynSys = simparams.dynSys;

if ~( strcmp(dynSys,'2bp') || strcmp(dynSys,'cr3bp') )
    error('Error: Dynamical system identifier text not recognized!')
end


if nargin == 3
    stm_initial = eye(6);
    stt_initial = zeros(6,6,6);
end

if nargin == 4
    stt_initial = zeros(6,6,6);
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
    end
    xstmstt_f = xstmstt_t(end,:);
    xfinal = reshape(xstmstt_f(1:6),6,1);
    stm = reshape(xstmstt_f(7:42),6,6);
    stt = reshape(xstmstt_f(43:end),6,6,6);
end






end