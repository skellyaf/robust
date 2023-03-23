% This program creates points on the circumference of a UNIT circle.
% It is the two dimensional analog of "sphere" the built-in Matlab
% command.  X and Y are arrays together containing n+1 x,y pairs 
% which describe a circle. See also "ellipse" which creates arrays
% of x,y pairs that describe an ellipse, and "ellipsoid" the three
% dimensional version of "ellipse".
%
% Syntax: [X,Y] = circle(n);


function [X,Y] = circle(n);
d = 2*pi/n;
a = 0:n;
a = d*a;
X = sin(a);
Y = cos(a);
