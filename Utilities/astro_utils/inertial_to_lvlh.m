function T = inertial_to_lvlh(x);

r = x(1:3,1);
v = x(4:6,1);

ir = r/norm(r);
h =  cross(r,v);
ih = h/norm(h);
iv = cross(ih,ir);
T = [iv ih ir]';
