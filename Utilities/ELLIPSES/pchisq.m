% Returns the probability associated with chi squared and n dimensions
% the probability is that a randomly chosen set of variables lie 
% within the value of chi squared.  (Note the CRC tables show the probab-
% ility of lying outside the value.) 
%
%  Input:  n - problem dimensionality
%          csq - chi squared
%
%  Output:  p - probability associated with input chi squared
%               and dimensionality
%
%  12/16/98  DPF  Borrowed from Dick Phillips
%                 Note that this function essentially solves the chi-square
%                 distribution (see p. 537, "CRC Standard Math Tables")

function p = chisq(n,csq);

if rem(n,2) > 0;
   f(1) = sqrt(pi/2) * erf(sqrt(csq/2)); h(1) = sqrt(pi/2); m = 3;
else
   f(2) = 1 - exp(-csq/2); h(2) = 1;    m = 4;
end;   

if n > 2;
   for i = m:2:n
      f(i) = (i-2)*f(i-2) - csq^((i-2)/2) * exp(-csq/2);
      h(i) = (i-2)*h(i-2);
   end;  
end;   

p = f(n) / h(n);
 


