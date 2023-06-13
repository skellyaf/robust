function [apse_constraint_eqn, apse_constraint_gradient] = apse_constraint_prev_vel_S(x_curr, x_prev_f, r_body, stm_prev, simparams)
%APSE_CONSTRAINT Constrains a state's position and velocity vectors to be
%orthogonal (ie, their dot product to be zero). The returned values are the
%constraint equation (dot product of the position and velocity at a node)
%and the gradient of how the 
%   Detailed explanation goes here


% Position of the spacecraft wrt system barycenter
r_sc = reshape(x_curr(1:3),3,1);
% Position of the spacecraft wrt the body we're looking for an apsis with
r_b_sc = r_sc - r_body;



% Velocity of the spacecraft in the rotating (synodic) frame at the end of
% the previous segment
v_sc_S = reshape(x_prev_f(4:6),3,1);
% Converting velocity to an inertial frame for the apsis calc


if strcmp(simparams.dynSys,'2bp')
    omega_SI = [0; 0; 0]; % Not a rotating frame

elseif strcmp(simparams.dynSys,'cr3bp')
    omega_SI = [0; 0; 1]; % Synodic frame rotation rate

end

% Inertial velocity
% v_sc_I = v_sc_S + cross(omega_SI, r_sc);
% Inertial velocity unit vector
% i_v_sc_I = v_sc_I / vecnorm(v_sc_I);

% Cross product matrix of the rotation rate
% omegaCross = crossMatrix(omega_SI);


% Dot product apse constraint equation
% apse_constraint_eqn = r_b_sc' * v_sc_I;
apse_constraint_eqn = r_b_sc' * v_sc_S;


%% Gradient equations

% wrt previous segment duration
dx_prev_f = stateDot(x_prev_f, simparams.mu, simparams.dynSys);
% apse_constraint_dtime_prev = r_b_sc' * (dx_prev_f(4:6) + cross(omega_SI,cross(omega_SI,r_sc)));
apse_constraint_dtime_prev = r_b_sc' * dx_prev_f(4:6);

% assembled wrt previous segment parameters
apse_constraint_gradient_prev = [r_b_sc'*stm_prev(4:6,:), apse_constraint_dtime_prev];
% r_b_sc'*stm_prev(:,4:6)'

IO = [eye(3,3), zeros(3,3)];

% wrt current segment
% apse_constraint_gradient_curr = [v_sc_I' * IO + r_b_sc'*omegaCross*IO, 0];
apse_constraint_gradient_curr = [v_sc_S' * IO, 0];


% assembled both prev and current
apse_constraint_gradient = [apse_constraint_gradient_prev, apse_constraint_gradient_curr];


end

