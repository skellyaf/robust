function [deltaV, deltaVs] = calcDeltaV(x,x_i_f,simparams)
%calcDeltaV Calculates impulsive delta V between segment final velocities
%and the initial velocities of the next segment
% maneuverSegments identify the segment number with a maneuver at the
% beginning of it




m = simparams.m; % the number of elements per segment
n = simparams.n; % the number of segments
maneuverSegments = simparams.maneuverSegments;

% Reshaping x so each column is a segment initial state and duration
x = reshape(x,m,n);

if ismember(0, maneuverSegments-1)
    vi = x(4:6,maneuverSegments);
    msm1 = maneuverSegments-1;
    vf = [simparams.x_init(4:6), x_i_f(4:6,msm1(2:end)) ];

elseif ismember(n+1, maneuverSegments)
    assert(0,'NEED TO ADD FUNCTIONALITY HERE')
else
    vi = x(4:6,maneuverSegments);
    vf = x_i_f(4:6,maneuverSegments-1);
    
end

deltaVs = vi-vf;
deltaV = sum(vecnorm(deltaVs));

end