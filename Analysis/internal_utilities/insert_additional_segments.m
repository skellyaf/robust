function [x_new, simparams] = insert_additional_segments(x, simparams, i, n)
% This function inserts n additional segments into segment i and returns a
% new x and simparams with updated parameters
% Logic is not incorporated to handle subdivision of the first or last
% segment, since those are coast segments, it shouldn't be required. But if
% other mods are made to the overall method, including some new logic for
% that will be required.


x = reshape(x,simparams.m,simparams.n);
nsv = simparams.m - 1;

% Seg i is long, divide it into 1+n segments (+n)

[~,x_i_t, t_i] = stateProp(x(1:nsv,i), x(simparams.m,i), simparams);
x_i_new = subdivide_segment(x_i_t, t_i, n+1);

x_new = zeros(simparams.m, size(x,2)+n);



x_new(:,1:i-1) = x(:,1:i-1);
x_new(:,i:i+n) = x_i_new;
x_new(:,i+n+1:end) = x(:,i+1:end);




% Re-do params

simparams.n = size(x_new,2);

simparams.maneuverSegments(simparams.maneuverSegments > i) = simparams.maneuverSegments(simparams.maneuverSegments > i) + n;
simparams.P_constrained_nodes = simparams.maneuverSegments(2:end);

if isfield(simparams,'tcm_nodes')

    simparams.tcm_nodes(simparams.tcm_nodes > i) = simparams.tcm_nodes(simparams.tcm_nodes > i) + n;
end

x_new = x_new(:);



end