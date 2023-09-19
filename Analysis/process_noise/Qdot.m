function [Qdot_out] = Qdot(Q, x, mu, Qt)


F = cr3bp_sFrame_nd_Amatrix(x, mu);

G = [zeros(3,3); eye(3,3)];

% stmdot = F*stm;

% Fab = cr3bp_stt2_tensor(x, mu);


% Qdot equations of motion
Qdot_out = F*Q + Q*F' + G*Qt*G';

% dFdx = einsum2(Fab, stm, [1, 2, 4], [4, 3], 3);

% Q2dot = tmult(F,Q2) +  tmult(Q2, F,[0 1]) + tmult(dFdx, Q) + tmult(Q, dFdx, [0 1]);


















end

