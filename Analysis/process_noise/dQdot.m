function [Q2dot] = dQdot(Q2, Q, x, stm, mu)


F = cr3bp_sFrame_nd_Amatrix(x, mu);

% stmdot = F*stm;

Fab = cr3bp_stt2_tensor(x, mu);


% Qdot equations of motion
% Qdot = F*Q + Q*F' + G*Qt*G';

% Q2dot = tmult(Fab, Q) + tmult(F,Q2) + tmult(Q, Fab, [0 1]) + tmult(Q2, F,[0 1]);


% pFab = ipermute(Fab, [1, 2, 4, 3]); %T
% pstm = ipermute(stm, [4, 3, 1, 2]); %A1
% dFdx = sum(  permute(  pFab.*pstm, [1, 2, 3, 4]  ), [4]   );
dFdx = einsum2(Fab, stm, [1, 2, 4], [4, 3], 3);
% dFdt = einsum2(Fab, stm, [1, 2, 4], [4, 3], 3);



Q2dot = tmult(F,Q2) +  tmult(Q2, F,[0 1]) + tmult(dFdx, Q) + tmult(Q, dFdx, [0 1]);


















end

