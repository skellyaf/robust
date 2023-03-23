function x = gen3segTraj(t1,t2,t3,simparams)

x0 = simparams.x_init;

% Propagate coast in first orbit from x0 to burn location
[x1, stm1] = stateStmProp(x0, t1, simparams);
r2_start = x1(1:3);
v2_start = x1(4:6);

% Propagate backwards coast in target orbit from x_target to arrival burn
[x3_start, stm3_T] = stateStmProp(simparams.x_target, -t3, simparams);
stm3 = stm3_T';
r3_start = x3_start(1:3);
v3_start = x3_start(4:6);


h_xfer = cross(r2_start, r3_start);
uh_xfer = h_xfer / norm(h_xfer);

if uh_xfer(3) < 0
    uh_xfer = -uh_xfer;
end

% v1 and v2 of xfer orbit using Lambert problem
[ v2_xfer_start, v2_xfer_end ] = lambert (simparams.mu*1e9 / 3600^2, t2*3600, r2_start*1000, r3_start*1000, uh_xfer);
v2_xfer_start = v2_xfer_start / 1e3 * 3600;
v2_xfer_end = v2_xfer_end / 1e3 * 3600;

x2_xfer_start = [r2_start; v2_xfer_start];

x = [x0; t1; x2_xfer_start; t2; x3_start; t3];


end