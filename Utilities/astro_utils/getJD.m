function JD=getJD(yr,mon,day,hr,min,s)

%JD=getJD(yr,mon,day,hr,min,s)
%Calculate the Julian date from a know calendar date
%input yr, mon,day,hr,min,sec

% Ref: Fundamentals of Astrodynamcis and Applications
%      D. A. Vallado, pg 67 (algorithm 2)

JD = 367*yr - fix(7*(yr+fix((mon+9)/12))/4) + fix(275*mon/9) ...
     + day + 1721013.5 + ((s/60+min)/60+hr)/24;
