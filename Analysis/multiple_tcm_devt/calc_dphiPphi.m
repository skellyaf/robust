function dphiPphi = calc_dphiPphi(phi, dphi, P, dP)

% depth_dT = size(dT,3);
% depth_dP = size(dP,3);
% assert(depth_dT == depth_dP);

dphiPphi = zeros(size(dphi));



t1 = tmult(dphi,P*phi');
t2 = tmult(phi,tmult(dP,phi'));
t3 = tmult(phi*P,dphi,[0 1]);

dphiPphi = t1 + t2 + t3;




end

