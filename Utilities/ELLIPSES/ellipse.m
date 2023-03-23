% function [X,Y] = ellipse(M,n);
%
% Program to create an array of points which describe the perimeter 
% an ellipse.  The input is a two dimensional covariance matrix.  
% The output is two arrays which together contain the x,y coordinates
% of the perimeter of the ellipse. The eigenvalues of the matrix are
% used to scale the x and y coordinates of the unit circle (obtained
% from the "circle" function).  The eigenvectors of the matrix are
% used to rotate the coordinate pairs to the proper orientation.
%
% Note to represent the locus equal probability density for a specific
% value of enclosed probability the points X and Y must be multiplied
% by the appropriate value of the sqrt of chisq.  E.g. for the ellipse
% which encloses 0.68 of the probability:
%
% fctr = sqrt(chisq(2,0.68));
% X    = fctr * X;
% Y    = fctr * Y;
% 
% Input: M the matrix representing the ellipse
%        n the number of points in the arrays X and Y
%

function [X,Y] = ellipse(M,n);
[ROT,EIGV] = eig(M);
[X,Y] = circle(n);

X = sqrt(EIGV(1,1)) * X; 
Y = sqrt(EIGV(2,2)) * Y; 

B = ROT * [X;Y];

X = B(1,:);
Y = B(2,:);
 
return

