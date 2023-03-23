function [xfinal, stm, stt, xstmstt_t] = stateStmSttProp_fixed(x_initial, t, t_step, mu, stm_initial, stt_initial)

% global options;


if nargin == 4
    stm_initial = eye(6);
    stt_initial = zeros(6,6,6);
end

if nargin == 5
    stt_initial = zeros(6,6,6);
end

negMult = 1;
if t < 0
    negMult = -1;
    t = -t;
end

t_span = 0:t_step:t;


if t_span(end) ~= t && t ~= 0
    t_span(end+1) = t;
end
t_span = t_span * negMult;




xstmstt_initial = [x_initial; stm_initial(:); stt_initial(:)];



if t ~= 0
    [xstmstt_t] = ode4(@r2bp_stt2_de, t_span, xstmstt_initial, mu);
    stt = reshape(xstmstt_t(end,43:end),6,6,6);    
    stm = reshape(xstmstt_t(end,7:42),6,6);    
    xfinal = reshape(xstmstt_t(end,1:6),6,1);
else
    stt = stt_initial;
    stm = stm_initial;
    xfinal = x_initial;
end







% [t, xstmstt_t] = ode45(@r2bp_stt2_de, [0,t], xstmstt_initial, options);
% xstmstt_f = xstmstt_t(end,:);
% xfinal = reshape(xstmstt_f(1:6),6,1);
% stm = reshape(xstmstt_f(7:42),6,6);
% stt = reshape(xstmstt_f(43:end),6,6,6);






end