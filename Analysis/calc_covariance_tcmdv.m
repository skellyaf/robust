function [P_target, tcm_dv_total, tcm_dv, P_i_minus, P_i_plus] = calc_covariance_tcmdv(x, t, t_s, stm_t, tcm_time, vel_disp_flag, P_i, simparams)
%calc_covariance_tcmdv Calculates the final dispersion covariance and the
%total delta V of the position corrections occuring along the nominal
%trajectory at each entry in the array tcm_time, as well as calculating the
%target dispersion covariance (P_target)


%     if sum(event_idx_logical)~=num_tcm
if length(tcm_time)~=length(unique(tcm_time))
    assert(0,'Error: you likely input duplicate TCM times in your tcm_time array');
end



m = simparams.m; % the number of elements per segment
n = simparams.n; % the number of segments

% Reshaping x so each column is a segment initial state and duration
x = reshape(x,m,n);

R_tcm = zeros(6,6);
R_tcm(4:6,4:6) = simparams.R; % TCM execution error covariance
G = [zeros(3,3); eye(3,3)];
R_dv = G*simparams.R_dv*G'; % nominal maneuver execution error covariance


%% Finding the nominal maneuvers, their times, and defining which events are TCMs

[event_times, event_is_tcm] = define_events(x(:), t, tcm_time, simparams);

% STM from beginning of trajectory to the end of the current stm history
% % This code is frequently called such that the end of the stm history is
% already the target. Additionally, the function define_events does not
% include the target as an event. Therefore, the final propagation outside
% the while loop brings the final "event" post covariance to the target.
stmN0 = stm_t(:,:,end);

% Empty structure to store the tcm_dv RSS magnitudes (not 3 sigma). It is
% one longer than the tcm_time because the velocity correction occurs at
% the end of the trajectory and isn't included in tcm_time.
num_tcm = length(tcm_time);
tcm_dv = zeros(1,length(tcm_time)+1);



%% Logic if no events (zero length)
% if length(tcm_time) == 0
if isempty(event_times)

    P_target = stmN0 * P_i * stmN0';
    P_i_minus = P_target;
    P_i_plus = nan;
    tcm_dv_total = 0;

else % Otherwise, if there are events, do:

    %% Setup for logic choices and STM portions
    % Test which segment has the correction
    
    event_idx_logical = logical(sum(t'==event_times', 1));    

    % Extract STMs for calculations
    % A tensor of STMs from the beginning of the trajectory to each correction
    stmC0 = stm_t(:,:,event_idx_logical);

    % Modified for events, not just TCMs

    % Plus one because a minus covariance for each event and the last is
    % going to be P at the target/end and there isn't an event for the
    % final / target time
    P_i_minus = zeros(6,6,length(unique(event_times))+1); 

    P_i_plus = zeros(6,6,length(unique(event_times))); % no plus one, it is the dispersion covariance after each event        
    
    j = 1;
    dup_cnt = 0;
    %% Compute the covariance only at the correction(s) and the final
%     for i = 1:length(event_times)
    while j <= length(unique(event_times))
        % j - for accessing variables that are being saved
        % i - for accessing event_times
        i = j + dup_cnt;
       
        % Propagate dispersion covariance from previous event to i

        if i == 1
            stmCiClast = stmC0(:,:,j);
        else
            stmCi0 = stmC0(:,:,j);
            stmCiClast = stmCi0 * stm0C;
        end 


        if j == 1 & event_idx_logical(1)
            P_i_minus(:,:,j) = P_i;
            stm0C = eye(6);
        else
            P_i_minus(:,:,j) = stmCiClast * P_i * stmCiClast';
            stm0C = invert_stm(stmC0(:,:,j), simparams);
            
        end


        
        stmNC = stmN0 * stm0C;

        % Is it a nominal maneuver, TCM, or both?
        % Check for both first
        if i == length(event_times)
            next_event_time = nan;
        else
            next_event_time = event_times(i+1);
        end
        if event_times(i) == next_event_time % if it is a nominal maneuver and TCM at the same time

            assert(xor(event_is_tcm(i), event_is_tcm(i+1)),'Some logic mishap happened...should only be here if one event is a TCM and one is a nominal manever, happening at the same time.');

            % Perform the TCM update   
        
            T = [-inv( stmNC(1:3,4:6) ) * stmNC(1:3,1:3), -eye(3)];
            N = [zeros(3,6); T];
            IN = eye(6) + N;
       
            P_i = IN * P_i_minus(:,:,j) * IN' + R_tcm + R_dv;
            P_i_plus(:,:,j) = P_i;
    
            P_tcm = T * P_i_minus(:,:,j) * T' + simparams.R; % with error in tcm_dv magnitude calc 
    %         P_tcm = T * P_i_minus(:,:,i) * T'; % without error in tcm_dv magnitude calc 

            tcm_idx = find(find(event_is_tcm)==i);
            if isempty(tcm_idx)
                tcm_idx = find(find(event_is_tcm)==i+1);
            end
            tcm_dv(tcm_idx) = sqrt(trace(P_tcm));


            dup_cnt = dup_cnt + 1;

%             i = i+1;

        elseif event_is_tcm(i) % otherwise, if it is a TCM
            % Perform the TCM update    
        
            T = [-inv( stmNC(1:3,4:6) ) * stmNC(1:3,1:3), -eye(3)];
            N = [zeros(3,6); T];
            IN = eye(6) + N;
       
            P_i = IN * P_i_minus(:,:,j) * IN' + R_tcm;
            P_i_plus(:,:,j) = P_i;
    
            P_tcm = T * P_i_minus(:,:,i) * T' + simparams.R; % with error in tcm_dv magnitude calc 
    %         P_tcm = T * P_i_minus(:,:,i) * T'; % without error in tcm_dv magnitude calc 
            tcm_idx = find(find(event_is_tcm)==i);
            tcm_dv(tcm_idx) = sqrt(trace(P_tcm));

        else % otherwise, it must be a nominal maneuver
            % Add maneuver execution error 
            P_i = P_i_minus(:,:,j) + R_dv;
            P_i_plus(:,:,j) = P_i;
            

        end

    
       

%         i = i+1;
        j = j+1;
        
    end
    
    % Propagate the rest of the way to the end of the trajectory
    % Need to do one more propagation because event_times does not

    
    P_target = stmNC * P_i * stmNC';
    P_i_minus(:,:,end) = P_target;


    if event_idx_logical(1)
        % No need to store a P_i_minus for the first "event" if it is the
        % initial covariance / nominal maneuver occurs at the initial time
        P_i_minus = P_i_minus(:,:,2:end);
    end



    if vel_disp_flag
        tcm_dv(end) = sqrt(trace(P_target(4:6,4:6)));
    else
        tcm_dv = tcm_dv(1:end-1);
    end

    tcm_dv_total = sum(tcm_dv);

end


end