function xstmdot = r2bp_stm_de(t, x, mu)

% global mu;

r = [x(1); x(2); x(3)];
v = [x(4); x(5); x(6)];
stm = reshape(x(7:42), 6,6);


f = [zeros(3,3), eye(3,3); 
    -mu/norm(r)^3 * eye(3,3), zeros(3,3)];

xdot = f * [r; v];

% Partial wrt state vector
A = r2bp_A_matrix([r; v], mu);

% STM propagation equation
stmdot = A*stm;

% Return
xstmdot = [xdot; reshape(stmdot, 36, 1)];
        

end