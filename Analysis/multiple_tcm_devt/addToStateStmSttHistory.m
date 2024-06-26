function [x_t, stm_t, stm_t_i, stt_t_i, t, t_s] = addToStateStmSttHistory(x_t, stm_t, stm_t_i, stt_t_i, t, t_s, add_times, simparams)
%addToStateStmSttHistory adds new time steps to the history data structures
%and numerically integrates the remainder of the data

add_times = sort(add_times);

% loop over each time to add
for i = 1:length(add_times)

    % only add the time if it already doesn't exist
    if ~ismember(add_times(i), t)
        

        % find all times before
        times_less_than_add_i = find(t < add_times(i));
        % find the index before
        idx_before_add_i = times_less_than_add_i(end);
        % find the segment of the index before  
        seg = t_s(idx_before_add_i);
        % Error checking - if the segment before and the segment after gap
        % a node
        assert(seg == t_s(idx_before_add_i + 1),'Need to incoroprate some additional logic if this fails');

        % state before
        x_before = x_t(idx_before_add_i,:)';
        % stm before
        stm_before = stm_t(:,:,idx_before_add_i);

        % find the corresponding stt_t_i{i} index
        stt_idx_before = find_stt_t_i_idx(t, t_s, stt_t_i, idx_before_add_i);
        % stt before the new time increment
        stt_before = stt_t_i{seg}(:,:,:,stt_idx_before);

        % calc delta_t
        delta_t = add_times(i) - t(idx_before_add_i);
        
        % propagate to the new state, stm, stt
%         [x_newTime, ~, ~, xstmstt_t_i, t_i] = stateStmSttProp(x_before, delta_t, simparams, stm_before, stt_before);
        [x_newTime, ~, ~, xstmstt_t_i, t_i] = stateStmSttProp(x_before, delta_t, simparams);
        
        % extract final time
        new_time = t_i(end);
        % extract final stm and stt
        xstmstt_new_time = xstmstt_t_i(end,:);
        stm_increment = reshape(xstmstt_new_time(7:42),6,6);
        stt_increment = reshape(xstmstt_new_time(43:end),6,6,6);
%         stm_new_time = reshape(xstmstt_new_time(7:42),6,6);
%         stt_new_time = reshape(xstmstt_new_time(43:end),6,6,6);
        

        % insert new state
        x_t = [x_t(1:idx_before_add_i,:); x_newTime'; x_t(idx_before_add_i+1:end,:)];

        % insert STM into time history
        stm_t_new = zeros(6,6,size(x_t,1));
        stm_t_new(:,:,1:idx_before_add_i) = stm_t(:,:,1:idx_before_add_i);
        stm_t_new(:,:,idx_before_add_i+1) = stm_increment * stm_before;
        stm_t_new(:,:,idx_before_add_i+2:end) = stm_t(:,:,idx_before_add_i+1:end);
        stm_t = stm_t_new;

        % insert STT into segment time history
        stt_t_i_new = zeros(6,6,6,size(stt_t_i{seg},4)+1);
        stt_t_i_new(:,:,:,1:stt_idx_before) = stt_t_i{seg}(:,:,:,1:stt_idx_before);
        stt_t_i_new(:,:,:,stt_idx_before+1) = tensorCombine(stm_before, stt_before, stm_increment, stt_increment);
        stt_t_i_new(:,:,:,stt_idx_before+2:end) = stt_t_i{seg}(:,:,:,stt_idx_before+1:end);
        stt_t_i{seg} = stt_t_i_new;

        stm_t_i_new = zeros(6,6,size(stm_t_i{seg},3)+1);
        stm_t_i_new(:,:,1:stt_idx_before) = stm_t_i{seg}(:,:,1:stt_idx_before);
        stm_t_i_new(:,:,stt_idx_before+1) = stm_increment * stm_t_i_new(:,:,stt_idx_before);
        stm_t_i_new(:,:,stt_idx_before+2:end) = stm_t_i{seg}(:,:,stt_idx_before+1:end);
        stm_t_i{seg} = stm_t_i_new;


        % insert into t history
        t = [t(1:idx_before_add_i); add_times(i); t(idx_before_add_i+1:end)];

        % insert into t_s history
        t_s = [t_s(1:idx_before_add_i); seg; t_s(idx_before_add_i+1:end)];


    


    end


end


end