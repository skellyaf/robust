function [xfinal, stm, xstm_t] = stateStmProp2t(xinitial, t1, t2,  stm_initial)

global options;

if nargin == 3
    stm_initial = eye(6);
end

xstm = [xinitial; reshape(stm_initial, 36, 1)];
if t2-t1 ~= 0
    [t, xstm_t] = ode45(@r2bp_stm_de, [t1,t2], xstm, options);
    stm = reshape(xstm_t(end,7:42),6,6);    
    xfinal = reshape(xstm_t(end,1:6),6,1);
else
    stm = stm_initial;
    xfinal = xinitial;
end
end