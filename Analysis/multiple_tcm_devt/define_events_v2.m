function [event_times, event_indicator] = define_events_v2(x, t, tcm_time, simparams)
%define_events Finds the nominal maneuvers, their times, and defining which events are TCMs

% In event_indicator:
%   0 = nominal maneuver
%   1 = TCM
%   2 = concurrent nominal maneuver and TCM

%%

m = simparams.m; % the number of elements per segment
n = simparams.n; % the number of segments

% Reshaping x so each column is a segment initial state and duration
x = reshape(x,m,n);



t_start = t(1);
% t_end = t(end);
t_end = max(t); % had a bug when the NLP tried to get the last seg to go backwards for some reason

maneuver_times = zeros(1,length(simparams.maneuverSegments));

for i = 1:length(simparams.maneuverSegments)
    maneuver_times(i) = sum(x(7,1:simparams.maneuverSegments(i)-1));
end

% Combine nominal DV and TCMs into one array of event times
event_times = [maneuver_times, tcm_time];
event_indicator = [zeros(1,length(maneuver_times)), ones(1,length(tcm_time))];

% Sort the event times and the indicator if it is a TCM
[event_times, sort_idx] = sort(event_times);
event_indicator = event_indicator(sort_idx);

% Eliminate the events before t_start and after t_end
keep_events = not(event_times < t_start | event_times > t_end); % inclusive at the beginning, exclusive at the end (we are already propagating to the end/target. include the maneuver execution error into the susbequent propagation) 
event_times = event_times(keep_events);
event_indicator = event_indicator(keep_events);


% If concurrent events, define as "2" in event_indicator and make it appear
% only once in event_times

[unique_event_times, unique_idx] = unique(event_times);

if length(unique_event_times) ~= length(event_times)

    are_unique = false(1,length(event_times));
    are_unique(unique_idx) = true;

    event_times = unique_event_times;


    event_indicator(find(~are_unique) - 1) = 2;
    event_indicator = event_indicator(unique_idx);
    


end


end