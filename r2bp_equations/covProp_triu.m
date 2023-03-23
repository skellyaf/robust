function [Pfinal, xfinal, xP_t] = covProp_triu(x_initial, Pinitial, t)
global mu options;

tri_u = [1, 7, 8, 13, 14, 15, 19, 20, 21, 22, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36]';
diag_el = [1, 8, 15, 22, 29, 36];

xP = [x_initial; Pinitial(tri_u)];

[~, xP_t] = ode45(@r2bp_covProp_triu, [0,t], xP, options);

Pfinal = zeros(6,6);
Pfinal(tri_u) = xP_t(end,7:end);

Pfinal = Pfinal' + Pfinal;
Pfinal(diag_el) = Pfinal(diag_el) / 2;

xfinal = reshape(xP_t(end,1:6),6,1);

% xP = [x_initial; reshape(Pinitial, 36, 1)];

% [t, xP_t] = ode45(@r2bp_covProp, [0,t], xP, options);
% Pfinal = reshape(xP_t(end,7:42),6,6);
% xfinal = reshape(xP_t(end,1:6),6,1);
end