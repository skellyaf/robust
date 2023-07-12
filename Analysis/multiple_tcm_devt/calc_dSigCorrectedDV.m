function dSig = calc_dSigCorrectedDV(i_DV, di_DV, Ptcm, dPtcm)

depth_di_DV = size(di_DV,2);
depth_dPtcm = size(dPtcm,3);
assert(depth_di_DV == depth_dPtcm);

% dSig = zeros(depth_di_DV,1);


sigDV = i_DV'*Ptcm*i_DV;

dTPT_t1 = tmult(di_DV',Ptcm*i_DV);
dTPT_t1 = dTPT_t1';
dTPT_t2 = tmult(i_DV',squeeze(tmult(dPtcm,i_DV)));
dTPT_t3 = tmult(i_DV'*Ptcm,di_DV,[0 0]);

% for i = 1:depth_di_DV
%     dSig(i) = trace( dTPT_t1(:,:,i) ) + trace( dTPT_t2(:,:,i) ) + trace( dTPT_t3(:,:,i) );
% end


dSig = dTPT_t1 + dTPT_t2 + dTPT_t3;

dSig = 1/2 * sigDV^(-1/2) * dSig;



end


% 
% 
% depth_dT = size(dT,3);
% depth_dPhi = size(dPhi,3);
% assert(depth_dT == depth_dPhi);
% 
% dvSigma_analytical_partial = zeros(depth_dT,1);
% 
% % dvSigma_analytical_partial_int_dT = tmult((T*phi*P*phi' + T*phi*P'*phi'), dT,[0 1]);
% dvSigma_analytical_partial_int_dT = tmult((T*phi*P*phi'), dT,[0 1]);
% 
% 
% 
% % dvSigma_analytical_partial_int_dPhi = tmult( ( T'*T*( phi*P + phi*P' ) ), dPhi, [0 1]);
% dvSigma_analytical_partial_int_dPhi = tmult( ( T'*T*phi*P  ), dPhi, [0 1]);
% 
% for i = 1:depth_dT
%     dvSigma_analytical_partial(i) = 2*( trace(dvSigma_analytical_partial_int_dT(:,:,i)) + trace(dvSigma_analytical_partial_int_dPhi(:,:,i)) );
% end