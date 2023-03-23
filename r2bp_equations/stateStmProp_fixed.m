function [xfinal, stm, xstm_t] = stateStmProp_fixed(x_initial, t, t_step, mu, stm_initial)

% global options;

if nargin == 4
    stm_initial = eye(6);
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

xstm = [x_initial; reshape(stm_initial, 36, 1)];
if t ~= 0
    [xstm_t] = ode4(@r2bp_stm_de, t_span, xstm, mu);
    stm = reshape(xstm_t(end,7:42),6,6);    
    xfinal = reshape(xstm_t(end,1:6),6,1);
else
    stm = stm_initial;
    xfinal = x_initial;
end
end