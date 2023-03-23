% function [X,Y,Z] = ellipsoid(M,n);
%
% This program creates an array of points which describe the surface
% of an ellipsoid.  The input is a three dimensional positive semi-
% definite (covariance for example) matrix.  The output is three
% two dimensional arrays which describe the x,y,z, coordinates of
% the surface of the ellipse.  Ellipsoid is the three-dimensional 
% version of ellipse.  It makes use of the Matlab "sphere" function
% the same way that "ellipse" makes use the "circle" function.  That
% is the eigenvalues are used to scale the x,y,z coordinates and the
% eigenvectors are used to rotate the ellipsoid to its proper 
% orientation.  The arrays X,Y,Z serve as inputs to the Matlab 3-D
% graphing functions mesh, surf, surfl, etc. 
%
% Note to represent the locus equal probability density for a specific
% value of enclosed probability the points X, Y, and Z must be multiplied
% by the appropriate value of the sqrt of chisq.  E.g. for the ellipsoid
% which encloses 0.68 of the probability:
%
% fctr = sqrt(chisq(3,0.68));
% X    = fctr * X;
% Y    = fctr * Y;
% Z    = fctr * Z;
%
% Inputs: 	M covariance matrix (3x3)
% 			n number of points used to describe ellipsoid
%

function [X,Y,Z] = ellipsoid(M,n);

[X,Y,Z] = sphere(n);
nps     = (n+1)^2;

[ROT,EIGV] = eig(M); 

X = sqrt(EIGV(1,1)) * X;
Y = sqrt(EIGV(2,2)) * Y;
Z = sqrt(EIGV(3,3)) * Z;

xrow = reshape(X,1,nps); 
yrow = reshape(Y,1,nps); 
zrow = reshape(Z,1,nps); 

B = [xrow;yrow;zrow]; 

B = ROT * B; 

X = reshape(B(1,:),n+1,n+1);
Y = reshape(B(2,:),n+1,n+1);
Z = reshape(B(3,:),n+1,n+1);

return
