function diDV = calc_diDVdDV(DV, dv_mag)


diDV = eye(3)/dv_mag - DV * DV' / dv_mag^3;






end

