function [stm_i, stt_i, x_i_f, x_t, stm_t, stt_t_i, t, t_s, stm_t_i] = createStateStmSttHistory(x, simparams)
%createStateStmHistory Returns the history of the state, STM, and STT via
%numerical integration of dynamics

% The output STT structure is different than the STM structure!!!
% stm_i is a 6x6xn tensor of STMs from the beginning to end of each n
% segments which extends from stm(tf,i , t0,i)
% x_i_f is a 6xn matrix of the final states at the end of each n segments
% x_t is a 6xj matrix of states along the trajectory that crosses all
% segments. j encompasses every automatically chosen numerical integration
% time step.
% stm_t is a 6x6xj tensor of STM histories, each 6x6 matrix index of which
% extends from stm(t(j), t0)
% stt_i is a 6x6x6xn STT structure that matches the structure of stm_i
% where each STT extends from the beginning to end of the i-th segment.
% stt_t_i is different than stm_t!!! stt_t_i is a 6x6x6xj {n} cell structure that
% contains the STT history for each individual segment, from the beginning
% of that segment. j represents the time history along an individual
% segment. n represents the number of segments. POTENTIALLY CONFUSING: it
% includes the 1:end element for each segment. the time vector does not so
% that it doesn't duplicate a time element. the time vector matches with
% 2:end of each STT segment history except the first STT cell which matches.

% t_s is an array that is j elements long that identifies which segment
% each time element occurs in. since the time at the beginning of each
% segment equals the time at the end of the previous segment, the times and
% histories are stored only once. for the cumulative state and STM history,
% the state and STM at the beginning of each segment is omitted along with
% the identical time element being omitted.

mu = simparams.mu;

% The correction is not part of the multi-segment trajectory vector x (as
% it has been in some instances)

m = simparams.m; % the number of elements per segment
n = simparams.n; % the number of segments

% Reshaping x so each column is a segment initial state and duration
x = reshape(x,m,n);

% Preallocate
stm_i = zeros(6,6,n);
stt_i = zeros(6,6,6,n);

% if nargin < 3
%     additional_t_steps = [];
%     include_more_t_steps = 0;
% else
%     additional_t_steps = sort(additional_t_steps);
%     include_more_t_steps = 1;
%     add_t_counter = 1;
% end

%% Propagate entire trajectory, save dynamics at each time step
% i = 1;
% run = 1;
% while run
for i = 1:n

%     % Get the total time along the traj at the beginning of the i-th seg 
%     if i == 0
%         t_start_total = 0;
%     else
%         t_start_total = sum(x(7,1:i-1));
%     end
% 
%     % if we want to include TCM times as integration steps
%     % and the TCM lies between the beginning and the end of the ith segment
%     if include_more_t_steps && t_start_total < additional_t_steps(add_t_counter) && additional_t_steps(add_t_counter) < t_start_total + x(7,i)
%         delta_t = additional_t_steps(add_t_counter)
%         add_t_counter = add_t_counter + 1;
%         add_to_i = 0;
% 
%     else
    
    x_i_initial = x(1:6,i);
    delta_t = x(7,i);
%         add_to_i = 1;
%     end


    %%%%%%%%%%%% WHAT IF DELTA_T IS NEGATIVE? HOW TO ALLOW THAT?  REQUIRES
    %%%%%%%%%%%% THOUGHT (one example resulted in convergence to an
    %%%%%%%%%%%% infeasible point)
    if delta_t ~= 0
%         [xfinal, stm, stt, xstmstt_t, ti]
        % Returns the stm and stt from the beginning to end of an
        % individual segment
        [x_i_final, stm_i(:,:,i), stt_i(:,:,:,i), xstmstt_t_i, t_i] = stateStmSttProp(x_i_initial, delta_t, simparams, eye(6), zeros(6,6,6));

        x_i_f(:,i) = x_i_final;

        if i == 1
            % State history
            x_t = xstmstt_t_i(:,1:6);
            % STM time history from the beginning of each segment to each
            % time index
            stm_t = reshape( xstmstt_t_i(:,7:42)',6,6,[] );
            % STT history - cell structure - from the beginning of each segment to each time
            stt_t_i{i} = reshape( xstmstt_t_i(:,43:258)',6,6,6,[] );
            % Mirroring the above with a STM history cell structure
            stm_t_i{i} = stm_t;
            % Time history
            t = t_i;
            % Segment corresponding to each time
            t_s = ones([length(t_i), 1]);


        else
            % Append to history structure
            % Time (exclude first time element, it duplicates the final of
            % the previous)
            t = [t; t_i(2:end) + sum(x(7,1:i-1))];

            % Corresponding segment
            t_s = [t_s; i*ones([length(t_i)-1, 1])];

            % State
            x_t = [x_t; xstmstt_t_i(2:end,1:6)];

            % Reshape new STM history tensor
            stm_t_i_curr = reshape( xstmstt_t_i(:,7:42)',6,6,[] );
            stm_t_i{i} = stm_t_i_curr;

            % Reshape new STT history tensor - into cell structure for
            % single segment STT history only (not from the beginning of
            % entire trajectory)
%             stt_t_i{i} = reshape( xstmstt_t_i(2:end,43:258)',6,6,6,[] );
            stt_t_i{i} = reshape( xstmstt_t_i(:,43:258)',6,6,6,[] );

            % Combine with previous final STT to continue history from beginning of trajectory

%             if isempty(stt_t)
%                 stt_t = stt_t_i;
%             else
% 
% %                 % The STM and STT from the beginning of traj to the end of
% %                 % last segment/beginning of current segment
% %                 stm1 = stm_t(:,:,end);
% %                 stt1 = stt_t(:,:,:,end);
% % 
% %                 for i = 2:size(stt_t_i,4)
% %                     stm2 = stm_t_i(:,:,i);
% %                     stt2 = stt_t_i(:,:,:,i);
% %                     stt_t(:,:,:,end+i-1) = tensorCombine(stm1, stt1, stm2, stt2);
% %                 end
%             end

            % Combine with previous final STM to continue history from beginning of trajectory            
            if isempty(stm_t)
                stm_t = stm_t_i_curr;
            else
                stm_t(:,:,end+1:end+size(stm_t_i_curr,3)-1) = tmult(stm_t_i_curr(:,:,2:end),stm_t(:,:,end));      
            end
            
            

        end

    else
        x_i_final = x_i_initial;
        x_i_f(:,i) = x_i_final;
        stm_i(:,:,i) = eye(6);
        stt_i(:,:,:,i) = zeros(6,6,6);
        stt_t_i{i} = zeros(6,6,6);
        stm_t_i{i} = eye(6);

        if i == 1
            
            stm_t = eye(6);
%             stt_t_i{i} = zeros(6,6,6);
            x_t = [x_i_initial'];
            t = 0;
            t_s = i;
        end
    end



%     i = i + add_to_i;
%     if i > size(x,2);
%         run = 0;
%     end

end







end