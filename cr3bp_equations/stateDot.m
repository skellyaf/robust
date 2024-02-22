function [xdot] = stateDot(x, simparams)

dynSys = simparams.dynSys;
mu = simparams.mu;

if strcmp(dynSys,'2bp')
    % Evaluation of 2 body dynamics
    xdot = r2bp_de(1, x, mu);
elseif strcmp(dynSys,'cr3bp')
    % Evaluation of 3 body dynamics
    xdot = cr3bp_sFrame_nd_de(1, x, mu);
elseif strcmp(dynSys,'br4bp_em')
    % Bicircular restricted four body problem, earth-moon frame
    assert(0,'not yet implemented in this function');

elseif strcmp(dynSys,'br4bp_sb1')
    % Bicircular restricted four body problem, Sun-B1 frame
    mub = simparams.mub;
    a4 = simparams.a4; 
    theta_em_dot = simparams.theta_em_dot;
    xdot = br4bp_sb1_state_de(1, x, mub, mu, a4, theta_em_dot);

end


end