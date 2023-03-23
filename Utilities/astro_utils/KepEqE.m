function E = KepEqE(M,e,tol)

if (-pi < M < 0) | (M > pi)
     Eold = M - e;
else
     Eold = M + e;
end

error = 1.0e12;

while error > tol

    Enew = Eold + (M-Eold + e*sin(Eold))/(1-e*cos(Eold));
    error = abs(Enew-Eold);
    Eold = Enew;

end

E = Enew;
