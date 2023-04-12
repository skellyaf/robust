function [event_times, event_is_tcm] = define_events(x, t, tcm_time, simparams)
%define_events Finds the nominal maneuvers, their times, and defining which events are TCMs

m = simparams.m; % the number of elements per segment
n = simparams.n; % the number of segments

% Reshaping x so each column is a segment initial state and duration
x = reshape(x,m,n);



t_start = t(1);
t_end = t(end);

maneuver_times = zeros(1,length(simparams.maneuverSegments));

for i = 1:length(simparams.maneuverSegments)
    maneuver_times(i) = sum(x(7,1:simparams.maneuverSegments(i)-1));
end

% Combine nominal DV and TCMs into one array of event times
event_times = [maneuver_times, tcm_time];
event_is_tcm = [false(1,length(maneuver_times)), true(1,length(tcm_time))];

% Sort the event times and the bool indicator if it is a TCM
[event_times, sort_idx] = sort(event_times);
event_is_tcm = event_is_tcm(sort_idx);

% Eliminate the events before t_start and after t_end
keep_events = not(event_times < t_start | event_times >= t_end); % inclusive at the beginning, exclusive at the end (we are already propagating to the end/target. include the maneuver execution error into the susbequent propagation) 
event_times = event_times(keep_events);
event_is_tcm = event_is_tcm(keep_events);



end