function [di_DVddt, di_DVdx, di_DVdDV] = calc_diDV(deltaVs, i, stm_i, x_i_f, target_leg, simparams)


DV = deltaVs(:,target_leg);

maneuverSegments = simparams.maneuverSegments;
nsv = simparams.nsv;

if ismember(i+1, maneuverSegments) || ismember(i, maneuverSegments)

    % i_DV    
    dv_mag = vecnorm(DV);
    di_DVdDV = eye(3)/dv_mag - DV * DV' / dv_mag^3;  
end


if i+1 == maneuverSegments(target_leg)
    dDVdx = - stm_i(4:6,:,i);
    
    xdot_xif_minus = stateDot(x_i_f(:,i), simparams);
    vdot_xif_minus = xdot_xif_minus(4:6);
    dDVdt = -vdot_xif_minus;
elseif i == maneuverSegments(target_leg)
    dDVdx = [zeros(3,3), eye(3,3), zeros(3,mod(nsv,6))];

    dDVdt = zeros(3,1); 
else
    dDVdt = zeros(3,1);
    dDVdx = zeros(3,nsv);
    di_DVdDV = zeros(3,3);
end


    
di_DVddt = di_DVdDV * dDVdt;
di_DVdx = di_DVdDV * dDVdx;



end

