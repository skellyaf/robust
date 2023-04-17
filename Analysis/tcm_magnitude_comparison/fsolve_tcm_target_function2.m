function [del_rf] = fsolve_tcm_target_function2(tcm_dv, xc, x_target, dt, simparams)



xr_c_guess = xc + [zeros(3,1); tcm_dv];

xr_final = stateProp(xr_c_guess, dt, simparams);

del_xf = xr_final - x_target;

del_rf = del_xf(1:3);


end