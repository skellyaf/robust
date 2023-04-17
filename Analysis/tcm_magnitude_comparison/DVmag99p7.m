function [dvmag99p7] = DVmag99p7(Cov)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Coefficients Poly44 (99.7 percent probability): (with 95% confidence bounds):
       p00 =       2.959; %  (2.955, 2.963)
       p10 =     0.04764; %  (0.01477, 0.08051)
       p01 =     0.04764; %  (0.01477, 0.08051)
       p20 =      0.3032; %  (0.193, 0.4135)
       p11 =     -0.4948; %  (-0.5799, -0.4097)
       p02 =      0.3032; %  (0.193, 0.4135)
       p30 =     -0.7522; %  (-0.9053, -0.5991)
       p21 =      0.6789; %  (0.5636, 0.7941)
       p12 =      0.6789; %  (0.5636, 0.7941)
       p03 =     -0.7522; %  (-0.9053, -0.5991)
       p40 =      0.8523; %  (0.7781, 0.9266)
       p31 =     -0.2829; %  (-0.3466, -0.2192)
       p22 =     -0.4157; %  (-0.4779, -0.3535)
       p13 =     -0.2829; %  (-0.3466, -0.2192)
       p04 =      0.8523; %  (0.7781, 0.9266)
       
%% Begin calculations 
sigs = sort(sqrt(eig(Cov)));
sigz = sigs(3);
x = sigs(1)/sigz;
y = sigs(2)/sigz;

%Linear model Poly44 (99.7 percent probability):
dvmag99p7 = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + p30*x^3 ... 
                    + p21*x^2*y + p12*x*y^2 + p03*y^3 + p40*x^4 + p31*x^3*y ...
                    + p22*x^2*y^2 + p13*x*y^3 + p04*y^4;
dvmag99p7 = sigz*dvmag99p7;
end

