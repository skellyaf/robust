function [tcm_gradient, tcm_gradient_r, tcm_gradient_v, dPCkminusdxi, dPCkminusddti] = calc_multiple_tcm_gradient_wQ(x, traj, tcm_time, tcm_idx, P_k_minus, dQ_k_km1_dxi, dQ_k_km1_ddti, deltaVs_nom, simparams)
%calc_tcm_gradient Computes and returns the analytical tcm gradient


%   Detailed explanation goes here

% Variable extraction
m = simparams.m;
n = simparams.n;
mu = simparams.mu;
% P_initial = simparams.P_initial;
R = simparams.R;
G = [zeros(3,3); eye(3,3)];
R_dv = G*simparams.R_dv*G';
% maneuverSegments = simparams.maneuverSegments;
num_r_corrections = length(tcm_idx);

% Reshaping x so each column is a segment initial state and duration
x = reshape(x,m,n);

% Velocity mapping matrix
Mv = G';

% traj
x_i_f = traj.x_i_f;
stm_i = traj.stm_i;
stm_t = traj.stm_t;
stm_t_i = traj.stm_t_i;
stt_t_i = traj.stt_t_i;
t = traj.t;
t_s = traj.t_s;

% Store a gradient vector for each individual TCM (r and v)
% Structure for each tcm gradient
tcm_gradient_r = zeros(n*m, num_r_corrections);
tcm_gradient_v = zeros(n*m, 1);


%% New construct - looping through "events" instead of just TCMs
% - An event is a TCM (1) or a node with a nominal maneuver (0), or a
% node where both occur (2).

[event_times, event_indicator] = define_events_v2(x(:), t, tcm_time, simparams);
event_idx_logical = logical(sum(traj.t'==event_times', 1));    
event_idxs = find(event_idx_logical);

% Find the first TCM target time
target_time = sum(x(7,1:simparams.P_constrained_nodes(1)-1));
target_node = simparams.P_constrained_nodes(1);
target_leg = 1;
previous_target_node = 1;

% event_time_not_tcm = event_times(event_indicator~=1);
% target_time = event_time_not_tcm(1);
target_idx = find(t==target_time);

num_events = length(event_times);


%% Nominal delta V unit vectors and magnitude
DV_norms = vecnorm(deltaVs_nom);
i_DVs = deltaVs_nom./DV_norms;


%% Loop through the corrections for each of the following sets of calcs

% Set up data structure for dPCkminusdxi: the partial derivative of the
% covariance matrix P (prior to correction C (the minus part)) with respect
% to each segment initial state
% For each correction there will be a 6x6x6 partial derivative of the
% covariance matrix P- with respect to each segment.
%%% to start - calculate for each correction and use each iteration, think
%%% there is no need to save all of them.

% Initialize partial derivative placeholders for each correction loop
dPCkminusdxi = zeros(6,6,6,length(event_indicator),n);
dPCkminusddti = zeros(6,6,length(event_indicator),n);
dTkdxi = zeros(3,6,6,length(event_indicator),n);
dTkddti = zeros(3,6,length(event_indicator),n);
% dPCkminusdxi = zeros(6,6,6,num_r_corrections,n);
% dPCkminusddti = zeros(6,6,num_r_corrections,n);
% dTkdxi = zeros(3,6,6,num_r_corrections,n);
% dTkddti = zeros(3,6,num_r_corrections,n);






% dLkdxi = zeros(3,6,6,num_r_corrections,n);
% dLkddti = zeros(3,6,num_r_corrections,n);
% dstmNCkdxi = zeros(6,6,6,num_r_corrections,n);
% dstmNCkddti = zeros(6,6,num_r_corrections,n);
% dPndxi = zeros(6,6,6,n);
% dPddti = zeros(6,6,n);


%%%%% TODO - FIX UP traj INPUT/OUTPUT
% traj.t = t;
% traj.t_s = t_s;
% traj.stm_t_i = stm_t_i;
% traj.Q_t_i = Q_t_i;
% traj.stt_t_i = stt_t_i;
% traj.dQ_t_i = dQ_t_i;
% traj.Q_i = Q_i;

% Process noise covariance at TCM k, accumulated from k-1 to k: Qbar_k_km1
% [Qbar_k_km1, dQ_k_km1_dxi] = calc_Q_events(traj, x, tcm_time, simparams);

% Derivative of process noise wrt i


tcm_idx = 0;

%% Construct gradients loop

% while k < num_events
% for k = 1:num_events - 1
for k = 1:num_events

    if event_indicator(k) > 0
        tcm_idx = tcm_idx + 1;
    end

    
    % Get the time the kth event occurs
%     tCk = tcm_time(k);
    tCk = event_times(k);

    % Get the index tCk occurs
%     tCk_idx = tcm_idx(k);
    tCk_idx = find(t==tCk);

    % Test which segments have a correction
%     tCk_seg = t_s(tCk_idx);
%     tk_seg = t_s(tCk_idx);
    
    % Get the nominal state at the correction
%     x_tcm = x_t(tCk_idx,:)';

    % If the first iteration, the "last" TCM occured at t=0
    if k == 1
        if simparams.start_P_growth_node > 1
            tClast = sum(x(7,1:simparams.start_P_growth_node-1));
            tClast_idx = find(t == tClast);
            
        else
        
            tClast = 0;
            tClast_idx = 1;
%             tlast = 0;
%             tlast_idx = 1;
        end
    end
    
    % Extract STMs for calculations    
    stmCkClast = dynCellCombine(t, t_s, tClast_idx, tCk_idx, simparams, stm_t_i);

   
    if event_indicator(k) > 0 % If a TCM occurs at k
        % Calculate the T matrix
        stmNCk = dynCellCombine(t, t_s, tCk_idx, target_idx, simparams, stm_t_i);
        Tk = [-inv( stmNCk(1:3,4:6) ) * stmNCk(1:3,1:3), -eye(3)];

        % Calculate the W and L matrices (used for TCMv calculation)
        Wk = [stmNCk(4:6,4:6), -stmNCk(4:6,1:3)];
        Lk = [Wk*Tk', zeros(3,3)];
    
        % Calculate N and IN matrices
        Nk = [zeros(3,6); Tk];
        INk = eye(6) + Nk;
    end   
    
    %% Loop through the segments / calculate partial derivatives / gradients
    

    for i = 1:n

        % FIND:
        
        % -- dTdxi
        % -- dWdxi
        % -- dLdxi
        % ---- The above 3 need the following, which vary by segment
        % ---- dstmC0dxi ((((replace with dstmCkClastdxi)))) --- requires new logic
        % ---- dstmFCdxi

        % For multiple corrections, need to also calculate the following:
        % -- dPkmdxi (the partial of the pre-kth-correction covariance wrt state i)
        % -----Needs: 
        % -----     dINdxi (includes T partial, which is already calculated)---T partial needs updating for multiple TCMs
        % -----     dstmCkCminusdxi (a complicated one, needs 5 logic options
        %               based on t0,i, tf,i, tCk, tClast)

        if i == 7 && k == 5
            debug=1;
        end

        
        [dstmCkClastdxi, dstmCkClastddti] = calc_dstmCkClast(x, x_i_f, t, t_s, stm_t, stm_i, stm_t_i, stt_t_i, tCk, tClast, i, simparams); % verified numerically


        if event_indicator(k) > 0
            % Calculate dstmNC and assemble dTk

            if isempty(tCk) || isempty(target_time)
                ppp=1;
            end


            [dstmNCkdxi, dstmNCkddti] = calc_dstmCkClast(x, x_i_f, t, t_s, stm_t, stm_i, stm_t_i, stt_t_i, target_time, tCk, i, simparams);

            dTkdxi(:,:,:,k,i) = T_partial(stmNCk, dstmNCkdxi);
            dTkddti(:,:,k,i) = T_partial(stmNCk, dstmNCkddti);


        end

        % Assemble dPCkminusdxi
        if k == 1 %%%% if the first correction
            % Propagate dP to first event from P_initial
% % % %             If P covariance growth doesn't start until the first nominal maneuver, don't grow dP 
% % % %             if traj.t_s(event_idxs(k)) == traj.t_s(event_idxs(k)+1) - 1 && traj.t_s(event_idxs(k)+1) == simparams.start_P_growth_node
% % % %                 Do nothing in this case, should be zeros
% % % %                 p=1;
% % % %                 dPCkminusdxi(:,:,:,k,i) = calc_dphiPphi(stmCkClast, dstmCkClastdxi, simparams.P_initial, zeros(6,6,6), zeros(6,6)) + dQ_k_km1_dxi(:,:,:,k,i);
% % % %                 dPCkminusddti(:,:,k,i) = calc_dphiPphi(stmCkClast, dstmCkClastddti, simparams.P_initial, zeros(6,6), zeros(6,6)) + dQ_k_km1_ddti(:,:,k,i);
% % % %             else
% % % %                 dPCkminusdxi(:,:,:,k,i) = calc_dphiPphi(stmCkClast, dstmCkClastdxi, simparams.P_initial, zeros(6,6,6), zeros(6,6)) + dQ_k_km1_dxi(:,:,:,k,i);
% % % %                 dPCkminusddti(:,:,k,i) = calc_dphiPphi(stmCkClast, dstmCkClastddti, simparams.P_initial, zeros(6,6), zeros(6,6)) + dQ_k_km1_ddti(:,:,k,i);
% % % %             end

            dPCkminusdxi(:,:,:,k,i) = calc_dphiPphi(stmCkClast, dstmCkClastdxi, simparams.P_initial, zeros(6,6,6), zeros(6,6)) + dQ_k_km1_dxi(:,:,:,k,i);
            dPCkminusddti(:,:,k,i) = calc_dphiPphi(stmCkClast, dstmCkClastddti, simparams.P_initial, zeros(6,6), zeros(6,6)) + dQ_k_km1_ddti(:,:,k,i);
        else

            dPClast_minus_dxi = dPCkminusdxi(:,:,:,k-1,i);
            dPClast_minus_ddti = dPCkminusddti(:,:,k-1,i);

            % Structure P_k_minus has the dispersion P matrix prior to each EVENT 
            PClast_minus = P_k_minus(:,:,k-1);           

            % Update differently if TCM vs. nominal maneuver event (or both):
            % The event that matters is the last event (number k-1) because the
            % update to the current event time either involves a TCM or doesn't
            % pending what the previous event was.

            if event_indicator(k-1) == 2
                % Add effect of both TCM and nominal maneuver (passing both
                % R matrices to the calc_dP function)
                % AND process noise

                dTlastdxi = dTkdxi(:,:,:,k-1,i);
                dTlastddti = dTkddti(:,:,k-1,i);

                dPCkminusdxi(:,:,:,k,i) = calc_dPckMinus(stmCkClast, dstmCkClastdxi, simparams.R + simparams.R_dv, Tlast, dTlastdxi, PClast_minus, dPClast_minus_dxi) + dQ_k_km1_dxi(:,:,:,k,i);
                dPCkminusddti(:,:,k,i) = calc_dPckMinus(stmCkClast, dstmCkClastddti, simparams.R + simparams.R_dv, Tlast, dTlastddti, PClast_minus, dPClast_minus_ddti) + dQ_k_km1_ddti(:,:,k,i);


            elseif event_indicator(k-1) == 1
                % TCM only event
                % AND process noise
                dTlastdxi = dTkdxi(:,:,:,k-1,i);
                dTlastddti = dTkddti(:,:,k-1,i);

                dPCkminusdxi(:,:,:,k,i) = calc_dPckMinus(stmCkClast, dstmCkClastdxi, simparams.R, Tlast, dTlastdxi, PClast_minus, dPClast_minus_dxi) + dQ_k_km1_dxi(:,:,:,k,i);
                dPCkminusddti(:,:,k,i) = calc_dPckMinus(stmCkClast, dstmCkClastddti, simparams.R, Tlast, dTlastddti, PClast_minus, dPClast_minus_ddti) + dQ_k_km1_ddti(:,:,k,i);
    
            elseif event_indicator(k-1) == 3
                % A corrected nominal maneuver
                % AND the effect of process noise

                dTlastdxi = dTkdxi(:,:,:,k-1,i);
                dTlastddti = dTkddti(:,:,k-1,i);

                dPCkminusdxi(:,:,:,k,i) = calc_dPckMinus(stmCkClast, dstmCkClastdxi, simparams.R_dv, Tlast, dTlastdxi, PClast_minus, dPClast_minus_dxi) + dQ_k_km1_dxi(:,:,:,k,i);
                dPCkminusddti(:,:,k,i) = calc_dPckMinus(stmCkClast, dstmCkClastddti, simparams.R_dv, Tlast, dTlastddti, PClast_minus, dPClast_minus_ddti) + dQ_k_km1_ddti(:,:,k,i);

            else

                % Incorporating the impact of nominal maneuver execution
                % error AND process noise
                dPCkminusdxi(:,:,:,k,i) = calc_dphiPphi(stmCkClast, dstmCkClastdxi, PClast_minus, dPClast_minus_dxi, R_dv) + dQ_k_km1_dxi(:,:,:,k,i);
                dPCkminusddti(:,:,k,i) = calc_dphiPphi(stmCkClast, dstmCkClastddti, PClast_minus, dPClast_minus_ddti, R_dv) + dQ_k_km1_ddti(:,:,k,i);    
            end            
        end       


        if event_indicator(k) == 1 || event_indicator(k) == 2
            % Assemble gradient w/TCM execution error incorporated in
            % sigma_tcm calc (simparams.R)
            % dxi
            tcm_gradient_r((i-1)*7+1:(i-1)*7+6,tcm_idx) = calc_dSigk(Tk, dTkdxi(:,:,:,k,i), P_k_minus(:,:,k), dPCkminusdxi(:,:,:,k,i), simparams.R);    
            % ddti
            tcm_gradient_r(i*7,tcm_idx) = calc_dSigk(Tk, dTkddti(:,:,k,i), P_k_minus(:,:,k), dPCkminusddti(:,:,k,i), simparams.R);
            

        elseif event_indicator(k) == 3
            if i < target_node && i >= previous_target_node

                % Ptcm
%                 Ptcm = Tk * P_k_minus(:,:,k) * Tk' + simparams.R;
                Ptcm = Tk * P_k_minus(:,:,k) * Tk';
                % dPtcm
                dPtcmdxi = calc_dABA(Tk, dTkdxi(:,:,:,k,i), P_k_minus(:,:,k), dPCkminusdxi(:,:,:,k,i));
                dPtcmddti = calc_dABA(Tk, dTkddti(:,:,k,i), P_k_minus(:,:,k), dPCkminusddti(:,:,k,i));

                i_DV = i_DVs(:,target_leg);
                [di_DVddt, di_DVdx] = calc_diDV(deltaVs_nom, i, stm_i, x_i_f, simparams.maneuverSegments, target_leg, simparams.dynSys, simparams.mu); 

                % Into the gradient structure
                % ddti_minus
                tcm_gradient_r((i)*7,tcm_idx) = calc_dSigCorrectedDV(i_DV, di_DVddt, Ptcm, dPtcmddti);
                % dxi_minus
                tcm_gradient_r((i-1)*7+1:(i-1)*7+6,tcm_idx) = calc_dSigCorrectedDV(i_DV, di_DVdx, Ptcm, dPtcmdxi);


    
                
            end



        end

        % The final tcmV gradient
        if k == num_events - 1
            % Target dispersion covariance
%             Pn = P_k_minus(:,:,end) + G * simparams.R * G'; % QUESTION: ADDING TCM EXECUTION ERROR DIRECTLY TO FINAL STATE...IS THAT RIGHT?
            Pn = P_k_minus(:,:,end); % QUESTION: ADDING TCM EXECUTION ERROR DIRECTLY TO FINAL STATE...IS THAT RIGHT?
            dPndxi = calc_dPckMinus(stmNCk, dstmNCkdxi, simparams.R, Tk, dTkdxi(:,:,:,k,i), P_k_minus(:,:,k), dPCkminusdxi(:,:,:,k,i));
            dPnddti = calc_dPckMinus(stmNCk, dstmNCkddti, simparams.R, Tk, dTkddti(:,:,k,i), P_k_minus(:,:,k), dPCkminusddti(:,:,k,i));   

            assert(event_indicator(k)==1,'Error: need to incorporate logic here for when the event prior to the final target is a nominal maneuver!');

            if simparams.correct_nominal_dvs
                % If the TSE vector addition reduction is applied
                i_DV = i_DVs(:,end); % Final DV
                [di_DVddt, di_DVdx] = calc_diDV(deltaVs_nom, i, stm_i, x_i_f, simparams.maneuverSegments, length(simparams.maneuverSegments), simparams.dynSys, simparams.mu); 


                % Into the gradient structure
                % ddti_minus
                tcm_gradient_v(i*7,1) = calc_dSigCorrectedDV(i_DV, di_DVddt, Pn(4:6,4:6), dPnddti(4:6,4:6));
                % dxi_minus
                tcm_gradient_v((i-1)*7+1:(i-1)*7+6,1) = calc_dSigCorrectedDV(i_DV, di_DVdx, Pn(4:6,4:6), dPndxi(4:6,4:6,:));




            else
                % Otherwise, if just directly cleaning up the remaining
                % velocity dispersion
    
                
    
    
                
          
    
                tcm_gradient_v((i-1)*7+1:(i-1)*7+6,1) = calc_dSigk(Mv, zeros(3,6,6), Pn, dPndxi, simparams.R);
                tcm_gradient_v(i*7,1) = calc_dSigk(Mv, zeros(3,6), Pn, dPnddti, simparams.R);


            end
           



        end

        
    end % end of for loop for each segment (inside the loop for each TCM k)


    tClast = tCk; % Updating tClast for the next iteration
    tClast_idx = tCk_idx;
%     tClast_seg = tCk_seg;
    


    
    if event_indicator(k) > 0
        % Updating Tlast, which is referenced the next iteration:
        Tlast = Tk;
    end

    % Assign new target time
    if k+1 <= num_events
    
        if event_times(k+1) == target_time
            future_targets = event_times(event_times>target_time & event_indicator~=1);
            target_time = min(future_targets);
            target_idx = find(t==target_time);    
    
            target_leg = target_leg + 1;
            if target_leg <= length(simparams.P_constrained_nodes)
                target_node = simparams.P_constrained_nodes(target_leg);
            end
        end
    end
        
%     if event_indicator(k+1) ~= 1
%     
%         if k < num_events
%             target_time = event_time_not_tcm(find(event_time_not_tcm==target_time)+1);
%             target_idx = find(t==target_time);     
%         end
%     
%     end

end % end of loop for each TCM


%% Calculate output
tcm_gradient = sum(tcm_gradient_r')' + tcm_gradient_v;


end