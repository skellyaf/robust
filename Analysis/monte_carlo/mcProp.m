function [x_t, stm_t, t, x_curr] = mcProp(x,event_times,event_dvs,simparams)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


stm_t = [eye(6)];
x_t = [x];
t = [0];

[~,sort_order] = sort(event_times);
event_times = event_times(sort_order);
event_dvs = event_dvs(:,sort_order);

x_curr = x;

for i = 1:length(event_times)
    if event_times(i) - t(end) > 0
        % Propagate to event at event time i
        [x_curr, stm_i, xstm_t_i, t_i] = stateStmProp(x_curr, event_times(i) - t(end), simparams);
        % add corresponding delta V
        dv = event_dvs(:,i);
        x_curr = x_curr + [zeros(3,1); dv];

        x_t = [x_t, xstm_t_i(2:end,1:6)'];
        t = [t; t_i(2:end) + t(end)];
        stm_t_i = reshape( xstm_t_i(:,7:42)',6,6,[] );
        stm_t(:,:,end+1:end+size(stm_t_i,3)-1) = tmult(stm_t_i(:,:,2:end),stm_t(:,:,end));

    else
        % Means we are already at event time i
        % Add corresponding delta v to x_curr
        dv = event_dvs(:,i);
        x_curr = x_curr + [zeros(3,1); dv];



    end

end

p=1;

end