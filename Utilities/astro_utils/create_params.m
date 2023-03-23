close all, clear all, clc;
% Script to create common orbital parameters and save as data file

G = 6.67430e-20; % (km^3/kg/s) Gravitational constant

% Parameters for the sun

sun.mu = 132712e6; % (km^3/s^2) Sun gravitational parameter
sun.mass = 1988500e24;% (kg) mass of the sun

%Parameters for Earth

earth.mu = 398600.4415; % (km^3/s^2) Earth gravitational parameter
earth.rad = 6378; % (km) equatorial radius of the earth
earth.a = 149.6e6; % (km) Semimajor axis of Earth's orbit around the sun
earth.mass = 5.9724e24; % (kg) mass of Earth
earth.rp = 147.092e6; % (km) radius of perihelion for Earth
earth.ra = 152.099e6; % (km) radius of apohelion for Earth

% Parameters for the Moon

moon.mu = 4900; % (km^3/s^2) Moon gravitational parameter
moon.a = 384400; % (km) Semimajor axis of moon orbit around Earth
moon.rp = 363300; % (km) Radius of perigee of the moon's orbit around Earth
moon.ra = 405500; % (km) Radius of apogee of the moons' orbit around Earth
moon.rad = 1737.4; % (km) volumetric mean radius of the Moon
moon.mass = .07346e24; % (kg) Mass of the Moon
moon.revolutionPeriod = 27.3214; % days
moon.synodicPeriod = 29.53; % days
moon.e = .0549; % eccentricity

% Parameters for Mars

mars.mu = 42828.38; % (km^3/s^2) Mars gravitational parameter
mars.a = 227.923e6; % (km) Semimajor axis of Mar's orbit around the sun
mars.rp = 206.617e6; % (km) Radius of perihelion
mars.ra = 249.229e6; % (km) Radius of apohelion
mars.rad = 3389.5; % (km) volumetric mean radius of Mars
mars.mass = .064171e24; % (kg) Mass of  Mars
mars.period = 689.980; % (days) Sidereal orbit period of Mars
mars.synodicPeriod = 779.94;
mars.e = .09341233; % Mars eccentricity
mars.i = 1.85061; % Mars inclination



%  Save to matlab data file called orbital_params.mat
save orbital_params.mat;
