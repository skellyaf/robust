function T=inertial2earthfixed(JD)

J2000 = 2451545.0; %Epoch of J2000

%Days past J2000 epoch
D = JD - J2000;

Wearth = 190.147 + 360.9856235*D;

w = Wearth*pi/180; % Convert to radians

T = [cos(w) sin(w) 0;
    -sin(w) cos(w) 0;
        0     0    1];