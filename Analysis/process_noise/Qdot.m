function [Qdot_out] = Qdot(Q, x, simparams)

dynSys = simparams.dynSys;
F = A_matrix(x, simparams);

if strcmp(dynSys,'2bp')
    G = [zeros(3,3); eye(3,3)];
elseif strcmp(dynSys,'cr3bp')
    G = [zeros(3,3); eye(3,3)];
elseif strcmp(dynSys,'br4bp_sb1')
    G = [zeros(3,3); eye(3,3); zeros(1,3)];
end


% stmdot = F*stm;

% Fab = cr3bp_stt2_tensor(x, mu);


% Qdot equations of motion
Qdot_out = F*Q + Q*F' + G*simparams.Qt*G';

% dFdx = einsum2(Fab, stm, [1, 2, 4], [4, 3], 3);

% Q2dot = tmult(F,Q2) +  tmult(Q2, F,[0 1]) + tmult(dFdx, Q) + tmult(Q, dFdx, [0 1]);


















end

