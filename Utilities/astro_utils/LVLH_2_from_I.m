function [T, Tinv] = LVLH_2_from_I(state,option);

%returns the inertial to LVLH transformation

r = state(1:3);
v = state(4:6);

if (size(r,2) == 3)
 r = r';
 v = v';
end

rmag = norm(r);
ir = r/rmag;
in = cross(ir,v);  in = in/norm(in);
iv = cross(in,ir); iv = iv/norm(iv);

T = [iv'; in'; ir'];
Tinv = T';

%Inertial to rotating LVLH frame
if (option == 1)
        B = T;
        Binv = Tinv;

	w = cross(r,v)/(r'*r);

        W = [ 0   -w(3)  w(2)
             w(3)   0   -w(1)
            -w(2)  w(1)   0  ];

        T =    [ B    zeros(3,3)
               -B*W,     B      ];

        Tinv = [ Binv    zeros(3,3)
                 W*Binv,   Binv      ];
else
%Inertial to quasi-inertial LVLH frame
     B = T;
     Binv = Tinv;
     T =    [ B         zeros(3,3)
             zeros(3,3)     B      ];

     Tinv = [ Binv       zeros(3,3)
              zeros(3,3)   Binv      ];

end
        



