function xdot = r2bp_de(t, x, mu)



r = [x(1); x(2); x(3)];
v = [x(4); x(5); x(6)];

f = [zeros(3,3), eye(3,3); 
    -mu/norm(r)^3 * eye(3,3), zeros(3,3)];

xdot = f * x;
        

end