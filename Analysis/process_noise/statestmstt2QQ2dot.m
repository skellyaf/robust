function [xDot] = statestmstt2QQ2dot(t, X, mu, Qt)

% unpack x
assert(length(X) == 510);

G = [zeros(3,3); eye(3,3)];

x = X(1);
y = X(2);
z = X(3);
vx = X(4);
vy = X(5);
vz = X(6);


stm = reshape(X(7:42),6,6);
stt2 = reshape(X(43:258),6,6,6);
Q = reshape(X(259:294),6,6);
% ADDED STT
Q2 = reshape(X(295:510),6,6,6);

% X dot equations of motion for CR3BP
r1 = sqrt(  (x + mu)^2 + y^2 + z^2  );
r2 = sqrt(  (x-1+mu)^2 + y^2 + z^2  );

ax = 2*vy + x - (1-mu)*(x+mu)/r1^3 - mu*(x-1+mu)/r2^3;
ay = -2*vx + y -(1-mu)*y/r1^3 - mu*y/r2^3;
az = -(1-mu)*z/r1^3 - mu*z/r2^3;


F = cr3bp_sFrame_nd_Amatrix(X(1:6), mu);

stmdot = F*stm;

Fab = cr3bp_stt2_tensor([x; y; z; vx; vy; vz], mu);


% Qdot equations of motion
Qdot = F*Q + Q*F' + G*Qt*G';

% Q2dot = tmult(Fab, Q) + tmult(F,Q2) + tmult(Q, Fab, [0 1]) + tmult(Q2, F,[0 1]);


pFab = ipermute(Fab, [1, 2, 4, 3]); %T
pstm = ipermute(stm, [4, 3, 1, 2]); %A1

dFdx = sum(  permute(  pFab.*pstm, [1, 2, 3, 4]  ), [4]   );





Q2dot = tmult(F,Q2) +  tmult(Q2, F,[0 1]) + tmult(dFdx, Q) + tmult(Q, dFdx, [0 1]);
% Q2dot = tmult(F, Q2) + tmult(tmult(Fab,stm),Q2) + tmult(Q, Fab, [0 1]) + tmult(Q2, tmult(Fab,stm),[0 1]);

stt2dot = tensorCombine(stm,stt2,F,Fab);
% stt_2_0 = tensorCombine(stm_1_0, stt_1_0, stm_2_1, stt_2_1);


% Q2dot = tmult(stt2dot, Q) + tmult(Fab, Q) + tmult(F,Q2) + tmult(Q, stt2dot, [0 1]) + tmult(Q2, F,[0 1]) + tmult(Q, Fab, [0 1]);


% Put it back together
xDot = [vx; vy; vz; ax; ay; az; stmdot(:); stt2dot(:); Qdot(:); Q2dot(:)];









% Q = reshape(x(7:42),6,6);




end

