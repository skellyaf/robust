function [di_DVddt, di_DVdx, di_DVdDV] = calc_diDVdDV(DV, i, stm_i, x_i_f, maneuverSegments, dynSys)


dv_mag = vecnorm(DV);
di_DVdDV = eye(3)/dv_mag - DV * DV' / dv_mag^3;


if ismember(i+1, maneuverSegments) || ismember(i, maneuverSegments)

    % i_DV


    
    dv_mag = DV_norms(target_leg);
    

    di_DVdDV = calc_diDVdDV(deltaVs_nom(:,target_leg), );     
end


if i+1 == simparams.maneuverSegments(target_leg)
    dDVdx = - stm_i(4:6,:,i);
    
    xdot_xif_minus = stateDot(x_i_f(:,i), mu, dynSys);
    vdot_xif_minus = xdot_xif_minus(4:6);
    dDVdt = -vdot_xif_minus;
elseif i == simparams.maneuverSegments(target_leg)
    dDVdx = [zeros(3,3), eye(3,3)];

    dDVdt = zeros(3,1); 
else
    dDVdt = zeros(3,1);
    dDVdx = zeros(3,6);
    di_DVdDV = zeros(3,3);
end


    
di_DVddt = di_DVdDV * dDVdt;



di_DVdx = di_DVdDV * dDVdx;



end

