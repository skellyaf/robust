function A = r2bp_A_matrix(x, mu)

r = [x(1); x(2); x(3)];
v = [x(4); x(5); x(6)];
norm_r = norm(r);
ur = r/norm_r;

A = zeros(6,6);
A(4:6,1:3) = -mu/norm_r^3 * (  eye(3) - 3*ur*ur'  );
A(1:3,4:6) = eye(3);

end