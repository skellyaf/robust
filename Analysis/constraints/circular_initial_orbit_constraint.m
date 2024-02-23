function [ceq, ceqGrad] = circular_initial_orbit_constraint(x_i_initial,simparams)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

x0 = simparams.x_init;
mu = simparams.mu; 
ceq = zeros(3,1);

if nargout > 1
    ceqGrad = zeros(3,7);
end

if strcmp(simparams.dynSys,'cr3bp')
    
    omega_SI_nd = [0; 0; 1];
    % Constrain initial radius
    r_S_fixed = x0(1:3);
    r_earth = [-mu; 0; 0];
    rmag_E_fixed = vecnorm(r_S_fixed - r_earth);
    
    v_S_fixed = x0(4:6);
    v_I_fixed = v_S_fixed + cross(omega_SI_nd, r_S_fixed);
    
    rmag_fixed = vecnorm(x0(1:3));
    vmag_fixed = vecnorm(x0(4:6));
    
    
    
    energy_fixed = vecnorm(v_I_fixed)^2 / 2 - simparams.mu_earth_nd / rmag_E_fixed;
    
    
    r_S_curr = x_i_initial(1:3);
    r_E_curr = r_S_curr - r_earth;
    rmag_S_curr = vecnorm(r_S_curr);
    rmag_E_curr = vecnorm(r_S_curr - r_earth);
    
    v_S_curr = x_i_initial(4:6);
    v_I_curr = v_S_curr + cross(omega_SI_nd, r_S_curr);
    i_v_I = v_I_curr / vecnorm(v_I_curr);
    
    vmag_S_curr = vecnorm(v_S_curr);
    
    energy_curr = vecnorm(v_I_curr)^2 / 2 - simparams.mu_earth_nd / rmag_E_curr;
    
    % Initial constraint 1: fixed radius from Earth
    ceq(1) = rmag_E_curr - rmag_E_fixed;
    
    % Initial constraint 2: fixed specific energy (2 body)
    ceq(2) = energy_curr - energy_fixed;
    
    % Initial constraint 3: orthogonal position and velocity vectors
    ceq(3) = r_E_curr' * v_I_curr;
    
    omegaCross = crossMatrix(omega_SI_nd);
    
    
    if nargout > 1
        i_r_E = reshape( (r_S_curr - r_earth) / rmag_E_curr , 3, 1);
    
        ceqGrad(1, :) = [i_r_E', zeros(1,3), 0];
        ceqGrad(2, :) = [v_I_curr' * omegaCross + simparams.mu_earth_nd / rmag_E_curr^2 * i_r_E', v_I_curr', 0];
        ceqGrad(3, :) = [v_I_curr' + r_E_curr'*omegaCross'*i_v_I, r_E_curr', 0];
    end

end

end