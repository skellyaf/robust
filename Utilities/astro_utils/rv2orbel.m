function [a,ec,i,W,w,nu] = rv2orbel(r,v,MU_BODY);

%Calculates classical orbital elements from r, v
%input: r,v, MU_BODY

% Ref: Fundamentals of Astrodynamcis and Applications
%      D. A. Vallado, pg 146 (algorithm 5)

rmag = norm(r);
v2 = v'*v;
K = [0 0 1]';

h = cross(r,v);
hmag = norm(h);

n = cross(K,h);
nmag = norm(n);

e = ((v2 - MU_BODY/rmag)*r-(r'*v)*v)/MU_BODY;
emag = norm(e);
ec = emag;
e2 = emag*emag;

E = v2/2-MU_BODY/rmag;
a = -MU_BODY/2/E;
p = a*(1-e2);
i = acos(h(3)/hmag);

W = acos(n(1)/nmag);
  if (n(2) < 0 )
        W = 2*pi - W;
  end
 
w = acos(n'*e/emag/nmag);
  if ( e(3) < 0 )
	w = 2*pi - w;
  end

nu = real(acos(e'*r/emag/rmag));
  if (r'*v < 0 )
	nu = 2*pi - nu;
  end

