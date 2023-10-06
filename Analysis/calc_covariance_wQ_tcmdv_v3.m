function [P_target, tcm_dv_total, tcm_dv, P_i_minus, P_i_plus] = calc_covariance_wQ_tcmdv_v3(x, traj, tcm_time, vel_disp_flag, deltaVs_nom, P_i, Q_k_km1, simparams)
%calc_covariance_tcmdv Calculates the final dispersion covariance and the
%total delta V of the position corrections occuring along the nominal
%trajectory at each entry in the array tcm_time, as well as calculating the
%target dispersion covariance (P_target)


%     if sum(event_idx_logical)~=num_tcm
if length(tcm_time)~=length(unique(tcm_time))
    assert(0,'Error: you likely input duplicate TCM times in your tcm_time array');
end

% For the comparison below, at the intersection of segments, a node is
% shared amongst both segments. The way it is currently set up has the end
% of the segment including the node. So the first node of a segment in traj.t_s
% will show up as a the last segment. This applies for all but the first
% segment.

m = simparams.m; % the number of elements per segment
n = simparams.n; % the number of segments

% Reshaping x so each column is a segment initial state and duration
x = reshape(x,m,n);

R_tcm = simparams.R;
R_dv = simparams.R_dv; % nominal maneuver execution error covariance
G = [zeros(3,3); eye(3,3)];


%% Finding the nominal maneuvers, their times, and defining which events are TCMs

[event_times, event_indicator] = define_events_v2(x(:), traj.t, tcm_time, simparams);

event_idx_logical = logical(sum(traj.t'==event_times', 1));    
event_idxs = find(event_idx_logical);

% STM from beginning of trajectory to the end of the current stm history
% % This code is frequently called such that the end of the stm history is
% already the target. Additionally, the function define_events does not
% include the target as an event. Therefore, the final propagation outside
% the while loop brings the final "event" post covariance to the target.
% stmN0 = stm_t(:,:,end);

% target_event_logic = event_indicator == 0 | event_indicator == 2 | event_indicator == 3;
% target_idxs = event_idxs(target_event_logic);

target_idx_counter = 1;
if simparams.P_constrained_nodes(1) == simparams.n + 1
    target_idx = length(traj.t);
else
    
    target_idx = find(traj.t_s==simparams.P_constrained_nodes(1),1) - 1; % minus one because the segment intersection is shared by both segments, but assigned to the previous segment in traj.t_s
    if isempty(target_idx)
        idxs_Pi_minus1 = find(traj.t_s == simparams.P_constrained_nodes(1)-1);
        target_idx = idxs_Pi_minus1(end);
    end
end

maneuver_segments = simparams.maneuverSegments;


% Empty structure to store the tcm_dv RSS magnitudes (not 3 sigma). It is
% one longer than the tcm_time because the velocity correction occurs at
% the end of the trajectory and isn't included in tcm_time.
% num_tcm = length(tcm_time);
% tcm_dv = zeros(1,length(tcm_time)+1);
tcm_dv = [];




%% Logic if no events (zero length)
% if length(tcm_time) == 0
if isempty(event_times)

%     idx_start_P_growth = 

    stmN0 = dynCellCombine(traj.t, traj.t_s, 1, target_idx, simparams, traj.stm_t_i);

    P_target = stmN0 * P_i * stmN0';
    P_i_minus = P_target;
    P_i_plus = nan;
    tcm_dv_total = 0;


else % Otherwise, if there are events, do:

    %% Setup for logic choices and STM portions
    % Test which segment has the correction
    
    

%     % Extract STMs for calculations
%     % A tensor of STMs from the beginning of the trajectory to each correction
%     stmC0 = stm_t(:,:,event_idx_logical);

    % Modified for events, not just TCMs

    % Plus one because a minus covariance for each event and the last is
    % going to be P at the target/end and there isn't an event for the
    % final / target time
    P_i_minus = zeros(6,6,length(unique(event_times))); 

    P_i_plus = zeros(6,6,length(unique(event_times)));        

    start_time = sum(x(7,1:simparams.start_P_growth_node-1));
    start_idx = find(traj.t == start_time);
%     start_node = simparams.start_P_growth_node;
    

    %% Compute the covariance only at the correction(s) and the final
    for i = 1:length(event_times)
        % Propagate dispersion covariance from previous event to i

%         if traj.t_s(event_idxs(i)+1) >= simparams.start_P_growth_node || simparams.start_P_growth_node == 1
            if i == 1
                idx_Clast = 1;
            else
                idx_Clast = event_idxs(i-1);            
            end 
    
            idx_Ci = event_idxs(i);   
    
            % If there is an event at the very beginning of the trajectory
            % OR we want the covariance growth to start at the first event
            if i == 1 && event_idx_logical(1) || event_idxs(i) == start_idx
%             if i == 1 && event_idx_logical(1) || traj.t_s(event_idxs(i)) == traj.t_s(event_idxs(i)+1) - 1 && traj.t_s(event_idxs(i)+1) == simparams.start_P_growth_node
                P_i_minus(:,:,i) = P_i;
            else
                Q_Ci = Q_k_km1(:,:,i);
                stmCiClast = dynCellCombine(traj.t, traj.t_s, idx_Clast, idx_Ci, simparams, traj.stm_t_i);
                P_i_minus(:,:,i) = stmCiClast * P_i * stmCiClast' + Q_Ci;            
            end
    
            if idx_Ci == target_idx && target_idx_counter < length(simparams.P_constrained_nodes)
    
                    target_idx_counter = target_idx_counter + 1;
                    if simparams.P_constrained_nodes(target_idx_counter) == simparams.n + 1
                        target_idx = length(traj.t);
                    else
                        target_idx = find(traj.t_s==simparams.P_constrained_nodes(target_idx_counter),1) - 1;
                    end
            end       
    
    
    
            
            
    
            % Is it a nominal maneuver (0), TCM (1), both simultaneously (2), or a corrected nominal maneuver (3)?
    
            if event_indicator(i) > 0 % If there is a TCM (by itself or combined)
    
                stmNC = dynCellCombine(traj.t, traj.t_s, idx_Ci, target_idx, simparams, traj.stm_t_i);
    
                % Perform the TCM update matrix calculations     
                T = [-inv( stmNC(1:3,4:6) ) * stmNC(1:3,1:3), -eye(3)];
                N = [zeros(3,6); T];
                IN = eye(6) + N;
    
                if event_indicator(i) == 3 % Corrected nominal maneuver
                    % Corrected nominal applies the nominal maneuver execution error
                    P_tcm = T * P_i_minus(:,:,i) * T'; % QUESTION!!!!!!! with or without error since the error is incorporated in the dispersion covariance update
                    P_i = IN * P_i_minus(:,:,i) * IN' + G*R_dv*G';
    
                    % Get i_dv
                    deltaV = deltaVs_nom(:, maneuver_segments == traj.t_s(event_idxs(i)) + 1 );
    
                    i_dv = deltaV / vecnorm(deltaV);
                    tcm_dv(end+1) = sqrt(i_dv' * P_tcm * i_dv);
                elseif event_indicator(i) == 1 % TCM by itself
                    P_tcm = T * P_i_minus(:,:,i) * T' + R_tcm; % with error in tcm_dv magnitude calc 
                    P_i = IN * P_i_minus(:,:,i) * IN' + G*R_tcm*G';
                    tcm_dv(end+1) = sqrt(trace(P_tcm));
                elseif event_indicator(i) == 2 % Concurrent nominal maneuver and TCM performed separately
                    P_tcm = T * P_i_minus(:,:,i) * T' + R_tcm; % with error in tcm_dv magnitude calc 
                    P_i = IN * P_i_minus(:,:,i) * IN' + G*R_tcm*G' + G*R_dv*G';
                    tcm_dv(end+1) = sqrt(trace(P_tcm));
                else
                    assert(0,'Something happened that we did not plan for!');
                end
    
                
                % Covariance update (TCM, nominal maneuver is turned on/off with rdv_mult)       
                P_i_plus(:,:,i) = P_i;
    
    
            else % otherwise, it must be a nominal maneuver only
                % Add maneuver execution error 
                P_i = P_i_minus(:,:,i) + G*R_dv*G';
                P_i_plus(:,:,i) = P_i;
    
    
%             end
        end

    
        
        
    end
    
    P_target = P_i_minus(:,:,end); 
    
    if vel_disp_flag
        tcm_dv(end+1) = sqrt(trace(P_target(4:6,4:6)));
    end

    

    tcm_dv_total = sum(tcm_dv);

end


end