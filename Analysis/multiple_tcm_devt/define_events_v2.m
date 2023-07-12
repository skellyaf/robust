function [event_times, event_indicator] = define_events_v2(x, t, tcm_time, simparams)
%define_events Finds the nominal maneuvers, their times, and defining which events are TCMs

% In event_indicator:
%   0 = nominal maneuver
%   1 = TCM
%   2 = concurrent nominal maneuver and TCM performed separately
%   3 = corrected nominal maneuver

%%

m = simparams.m; % the number of elements per segment
n = simparams.n; % the number of segments

% Reshaping x so each column is a segment initial state and duration
x = reshape(x,m,n);


% if nargin == 5
%     t_start = t(range(1));
%     t_end = t(range(2));
% else
%     t_start = t(1);
%     % t_end = t(end);
%     t_end = max(t); % had a bug when the NLP tried to get the last seg to go backwards for some reason
% end

t_start = t(1);
% t_end = t(end);
t_end = max(t); % had a bug when the NLP tried to get the last seg to go backwards for some reason

maneuver_times = zeros(1,length(simparams.maneuverSegments));

for i = 1:length(simparams.maneuverSegments)
    maneuver_times(i) = sum(x(7,1:simparams.maneuverSegments(i)-1));
end

% Combine nominal DV and TCMs into one array of event times
event_times = [maneuver_times, tcm_time];

% Determine maneuver event indicator:
if simparams.correct_nominal_dvs
    maneuver_event_indicator = [3*ones(1,length(maneuver_times)-1), 0];
else
    maneuver_event_indicator = zeros(1,length(maneuver_times));
end

% TCM event indicator:
tcm_event_indicator = ones(1,length(tcm_time));

% Create combined event indicator:
event_indicator = [maneuver_event_indicator, tcm_event_indicator];

% Sort the event times and the indicator if it is a TCM
[event_times, sort_idx] = sort(event_times);
event_indicator = event_indicator(sort_idx);

% Eliminate the events before t_start and after t_end
keep_events = not(event_times < t_start | event_times >= t_end); % inclusive at the beginning, exclusive at the end (we are already propagating to the end/target. include the maneuver execution error into the susbequent propagation) 
event_times = event_times(keep_events);
event_indicator = event_indicator(keep_events);


% If concurrent events, define as "2" in event_indicator and make it appear
% only once in event_times

%% Handing duplicate events
% If a 0 (nominal maneuver) and a 1 (TCM) happen concurrently, combine them
% to be a 2 (concurrent nominal and TCM)

% However, if a corrected nominal maneuver (3) and a TCM (1) occur
% simultaneously, ....... %%%%TBD%%%%% PLACEHOLDER / ASSERT FOR FUTURE WORK
% 

% for testing, comment out later:
% event_times     = [0  0 .2 .3 .4 .4 .5]
% event_indicator = [3  1  1  0  1  0  0]

[unique_event_times, unique_idx, ic] = unique(event_times);

if length(unique_event_times) ~= length(event_times)


    are_unique = false(1,length(event_times));
    are_unique(unique_idx) = true;

    event_times = unique_event_times;

    not_unique = ~are_unique;


    tst_condition = event_indicator(not_unique) + event_indicator(find(not_unique)-1);

    condition_1 = tst_condition==1; % nominal maneuver and TCM
    condition_2 = tst_condition==4; % corrected nominal maneuver and TCM

    if sum(condition_2) > 0
        assert(0,'Need logic here for what happens when a corrected nominal maneuver is combined with a TCM!');
    end

    apply_condition = condition_1 * 2 + condition_2 * 3;

    event_indicator(find(not_unique) - 1) = apply_condition;



    event_indicator = event_indicator(unique_idx);
    


end


end