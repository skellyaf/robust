function [A] = A_matrix(x, simparams)


mu = simparams.mu;
dynSys = simparams.dynSys;


if strcmp(dynSys,'2bp')
    % Evaluation of 2 body dynamics
    A = r2bp_A_matrix(x, mu);
elseif strcmp(dynSys,'cr3bp')
    % Evaluation of 3 body dynamics
    A = cr3bp_sFrame_nd_Amatrix(x, mu);
elseif strcmp(dynSys,'br4bp_sb1')
    A = br4bp_sb1_Amatrix(x, simparams.mub, mu, simparams.a4);
end



end