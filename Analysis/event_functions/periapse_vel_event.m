function [rdotv,isterminal,direction] = periapse_vel_event(t,x,mu)
%%%%%%%%%%%

r_earth = [-mu, 0, 0]';

r_e2x = x(1:3) - r_earth;

vel = x(4:6);

rdotv = r_e2x'*vel;



% v = vecnorm(x(4:6));
% if t < -2
%     r_earth = [1, 0, 0]';
%     r_e2x = x(1:3) - r_earth;
%     
%     mag_e2x = norm(r_e2x);
%     
%     zero_val = -r_e2x'/mag_e2x * x(4:6);
% else 
%     zero_val = 1;
% end

% 
% X_e2man(i,:) = X_manifold_pos(i,:) - r_earth;
% X_e2man_mag(i) = norm(X_e2man(i,:));
% 
% X_e2_derivative(i) = - X_e2man(i,:)/(X_e2man_mag(i)) * X_manifold(i,4:6)';



isterminal = 1;  
direction = 1;