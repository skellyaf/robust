function K = KepEqK(l, P1, P2, tol)


Kold = l ;

error = 1.0e12;

i=0;
while error > 100*tol

    Knew = Kold + (l-Kold - P1*cos(Kold) + P2*sin(Kold))/(1-P1*sin(Kold)-P2*cos(Kold));
    error = abs(Knew-Kold);
    Kold = Knew;

i=i+1;
if i > 500
    keyboard
end

end

K = Knew;