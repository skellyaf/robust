function dSig = calc_dSigkdxi(T, dT, P, dP)

depth_dT = size(dT,3);
depth_dP = size(dP,3);
assert(depth_dT == depth_dP);

dSig = zeros(depth_dT,1);


sigk = sqrt(trace(T*P*T'));

dTPT_t1 = tmult(dT,P*T');
dTPT_t2 = tmult(T,tmult(dP,T'));
dTPT_t3 = tmult(T*P,dT,[0 1]);

for i = 1:depth_dT
    dSig(i) = trace( dTPT_t1(:,:,i) ) + trace( dTPT_t2(:,:,i) ) + trace( dTPT_t3(:,:,i) );
end

dSig = 1/2 * sigk^-1 * dSig;



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