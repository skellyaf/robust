function dABA = calc_dABA(A, dA, B, dB)




t1 = tmult(dA,B*A');
t2 = tmult(A,tmult(dB,A'));
t3 = tmult(A*B,dA,[0 1]);




dABA = t1 + t2 + t3;




end

