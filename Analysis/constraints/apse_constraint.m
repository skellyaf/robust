function [apse_constraint_eqn, apse_constraint_gradient] = apse_constraint(x_curr, r_body)
%APSE_CONSTRAINT Constrains a state's position and velocity vectors to be
%orthogonal (ie, their dot product to be zero). The returned values are the
%constraint equation (dot product of the position and velocity at a node)
%and the gradient of how the 
%   Detailed explanation goes here


% Position of the spacecraft wrt system barycenter
r_sc = reshape(x_curr(1:3),3,1);
% Position of the spacecraft wrt the body we're looking for an apsis with
r_b_sc = r_sc - r_body;



% Velocity of the spacectraft in the rotating (synodic) frame
v_sc_S = reshape(x_curr(4:6),3,1);
% Converting velocity to an inertial frame for the apsis calc


%%%% TODO: ADD AN IF STATEMENT FOR DYNAMICAL SYSTEM...THE BELOW ONLY
%%%% APPLIES TO CR3BP
omega_SI = [0; 0; 1]; % Synodic frame rotation rate

% Inertial velocity
v_sc_I = v_sc_S + cross(omega_SI, r_sc);
% Inertial velocity unit vector
i_v_sc_I = v_sc_I / vecnorm(v_sc_I);

% Cross product matrix of the rotation rate
omegaCross = crossMatrix(omega_SI);


% Dot product apse constraint equation
apse_constraint_eqn = r_b_sc' * v_sc_I;


% Gradient equation
apse_constraint_gradient = [v_sc_I' + r_b_sc'*omegaCross, r_b_sc', 0];



end

