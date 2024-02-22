function [P_target, tcm_dv_total, tcm_dv, P_i_minus, P_i_plus] = calc_covariance_wQ_tcmdv(x, traj, tcm_time, vel_disp_flag, deltaV, P_i, simparams)
%calc_covariance_tcmdv Calculates the final dispersion covariance and the
%total delta V of the position corrections occuring along the nominal
%trajectory at each entry in the array tcm_time, as well as calculating the
%target dispersion covariance (P_target)


% traj.x_t, t, traj.t_s, traj.stm_t

%     if sum(event_idx_logical)~=num_tcm
if length(tcm_time)~=length(unique(tcm_time))
    assert(0,'Error: you likely input duplicate TCM times in your tcm_time array');
end

% Make sure only a single deltaV vector was passed
assert(size(deltaV,2)==1 || size(deltaV,2)==2);

m = simparams.m; % the number of elements per segment
n = simparams.n; % the number of segments
nsv = simparams.nsv;


% Reshaping x so each column is a segment initial state and duration
x = reshape(x,m,n);

R_tcm = simparams.R;
R_dv = simparams.R_dv; % nominal maneuver execution error covariance
G = [zeros(3,3); eye(3,3); zeros(mod(nsv,6),3)];


%% Finding the nominal maneuvers, their times, and defining which events are TCMs

% [event_times, event_is_tcm] = define_events(x(:), t, tcm_time, simparams);

[event_times, event_indicator] = define_events_v2(x(:), traj.t, tcm_time, simparams);

% STM from beginning of trajectory to the end of the current stm history
% % This code is frequently called such that the end of the stm history is
% already the target. Additionally, the function define_events does not
% include the target as an event. Therefore, the final propagation outside
% the while loop brings the final "event" post covariance to the target.
target_idx = size(traj.stm_t,3);
% stmN0 = traj.stm_t(:,:,end);

% stmN0 = dynCellCombine(t, traj.t_s, start_idx, target_idx, simparams, traj.stm_t_i);


% Empty structure to store the tcm_dv RSS magnitudes (not 3 sigma). It is
% one longer than the tcm_time because the velocity correction occurs at
% the end of the trajectory and isn't included in tcm_time.
% tcm_dv = zeros(1,length(tcm_time)+1);
tcm_dv = [];




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
    
    event_idx_logical = logical(sum(traj.t'==event_times', 1));    
    event_idxs = find(event_idx_logical);

    % Extract STMs for calculations
    % A tensor of STMs from the beginning of the trajectory to each correction
%     stmC0 = traj.stm_t(:,:,event_idx_logical);


    % Modified for events, not just TCMs

    % Plus one because a minus covariance for each event and the last is
    % going to be P at the target/end and there isn't an event for the
    % final / target time
    P_i_minus = zeros(nsv,nsv,length(unique(event_times))+1); % Plus one because it captures the end of the trajectory portion, which isn't defined by an event

    P_i_plus = zeros(nsv,nsv,length(unique(event_times))); % no plus one, it is the dispersion covariance after each event    


%     Q_Clast = zeros(6,6);
    

    %% Compute the covariance only at the correction(s) and the final
    for i = 1:length(event_times)
       
        % Propagate dispersion covariance from previous event to i

        if i == 1
%             stmCiClast = stmC0(:,:,i);
            
            idx_Ci = event_idxs(i);
            idx_Clast = 1;
            stmCiClast = dynCellCombine(traj.t, traj.t_s, idx_Clast, idx_Ci, simparams, traj.stm_t_i);
            Q_Ci = traj.Q_t(:,:,idx_Ci);


            
        else
            idx_Ci = event_idxs(i);
            idx_Clast = event_idxs(i-1);
%             stmCi0 = stmC0(:,:,i);
%             stmCiClast = stmCi0 * stm0C;
            stmCiClast = dynCellCombine(traj.t, traj.t_s, idx_Clast, idx_Ci, simparams, traj.stm_t_i);
            Q_Ci = traj.Q_t(:,:,idx_Ci) - stmCiClast * traj.Q_t(:,:,idx_Clast) * stmCiClast';
        end 



        if i == 1 && event_idx_logical(1)
%         if i == 1 && event_idx_logical(1) || traj.t_s(event_idxs(i)) == traj.t_s(event_idxs(i)+1) - 1 && traj.t_s(event_idxs(i)+1) == simparams.start_P_growth_node
            P_i_minus(:,:,i) = P_i;
%             stm0C = eye(nsv);
        else
%             x_Clast = traj.x_t(idx_Clast,:)';
%             dt_CiClast = traj.t(idx_Ci) - traj.t(idx_Clast);
%             Q_Clast = zeros(6,6);



            %%%%%%% PROCESS NOISE %%%%%%%%
            % Want the process noise generated since the last time it was
            % included in the state covariance (which resets the process
            % noise to a 6x6 matrix of zeros) or from the beginning of the
            % trajectory. 

            % Have: idx_Ci, the index of the current event (where we want the
            % process noise effect to end); idx_Clast (where we want the
            % process noise effect to start).


%             Q_Ci = Qcombine(traj, idx_Clast, idx_Ci, simparams.dynSys);


            % Qcombine_2 update!
            % traj now includes Qbar_events and Qbart, each of which are
            % the Qbar at specific times along the trajectory that
            % correspond to growth from the previous Qbar zeroing event.

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%             [~, Q_Ci2] = stateQProp(x_Clast, dt_CiClast, simparams, Q_Clast);

%             Q_Ci = traj.Qbart(:,:,idx_Ci);
            
            P_i_minus(:,:,i) = stmCiClast * P_i * stmCiClast' + Q_Ci;
%             stm0C = invert_stm(stmC0(:,:,i), simparams);
        end


        
%         stmNC = stmN0 * stm0C;
        stmNC = dynCellCombine(traj.t, traj.t_s, idx_Ci, target_idx, simparams, traj.stm_t_i);

        % Is it a nominal maneuver (0), TCM (1), both at the same time (2), or a corrected nominal DV (3)?

        if event_indicator(i) > 0 % If there is a TCM (by itself or combined)

            % Perform the TCM update matrix calculations     
            T = [-inv( stmNC(1:3,4:6) ) * stmNC(1:3,1:3), -eye(3), zeros(3,mod(nsv,6))];
            N = [zeros(3,nsv); T; zeros(mod(nsv,6),nsv)];
            IN = eye(nsv) + N;

            if event_indicator(i) == 3 % Corrected nominal maneuver
                assert(simparams.correct_nominal_dvs); % Should only be here if this flag is on
                % Corrected nominal applies the nominal maneuver execution error
%                 P_tcm = T * P_i_minus(:,:,i) * T' + R_tcm; %%%%% QUESTION: APPLYING TCM R TO THE EXECUTION COST; APPLYING DV R TO THE DISPERSION...WHATS RIGHT? 
                P_tcm = T * P_i_minus(:,:,i) * T'; %%%%% QUESTION: APPLYING TCM R TO THE EXECUTION COST; APPLYING DV R TO THE DISPERSION...WHATS RIGHT? 
                P_i = IN * P_i_minus(:,:,i) * IN' + G*R_dv*G';

                % Get i_dv
%                 deltaV = deltaVs_nom(:, maneuver_segments == traj.t_s(event_idxs(i)) + 1 );

                i_dv = deltaV(:,1) / vecnorm(deltaV(:,1));
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

            % Covariance update    
            
            P_i_plus(:,:,i) = P_i;
  

        else % otherwise, it must be a nominal maneuver only

            
            % Add maneuver execution error 

            P_i = P_i_minus(:,:,i) + G*R_dv*G';
            P_i_plus(:,:,i) = P_i;

%             Q_Clast = Q_Ci;
            

        end

    
       

        
    end
    
    % Propagate the rest of the way to the end of the trajectory
    % Need to do one more propagation because event_times does not contain
    % the target

%     x_Ci = traj.x_t(idx_Ci,:)';
%     dt_targ_Ci = traj.t(end) - traj.t(idx_Ci);
%     Q_Ci = zeros(6,6);
% 
%     [~, Q_targ] = stateQProp(x_Ci, dt_targ_Ci, simparams, Q_Ci);
    

    Q_targ = traj.Q_t(:,:,end) - stmNC * traj.Q_t(:,:,idx_Ci) * stmNC';
    
    P_target = stmNC * P_i * stmNC' + Q_targ;
    P_i_minus(:,:,i+1) = P_target;


%     if event_idx_logical(1)
%         % No need to store a P_i_minus for the first "event" if it is the
%         % initial covariance / nominal maneuver occurs at the initial time
%         P_i_minus = P_i_minus(:,:,2:end);
%     end






    if vel_disp_flag
        if simparams.correct_nominal_dvs
            % Use the first order TSE expansion savings, if flagged
%             [~, DVs_nom] = calcDeltaV(x, traj.x_i_f, traj.stm_i, simparams)
%             DV_final = DVs_nom(:,end);
            i_dvf = deltaV(:,2) / vecnorm(deltaV(:,2));
            
            tcm_dv(end+1) = sqrt(i_dvf' * (P_target(4:6,4:6)) * i_dvf);
        else
            % If not flagged to perform simultaneous correction, just clean
            % up the remaining velocity dispersion
            tcm_dv(end+1) = sqrt(trace(P_target(4:6,4:6) + R_tcm));
        end
    end

    tcm_dv_total = sum(tcm_dv);

end


end