function [xdot] = stateDot(x, mu, dynSys)

if strcmp(dynSys,'2bp')
    % Evaluation of 2 body dynamics
    xdot = r2bp_de(1, x, mu);
elseif strcmp(dynSys,'cr3bp')
    % Evaluation of 3 body dynamics
    xdot = cr3bp_sFrame_nd_de(1, x, mu);
end


end