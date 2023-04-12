function [tcm_gradient, tcm_gradient_r, tcm_gradient_v] = calc_multiple_tcm_gradient(x, x_t, x_i_f, stm_i, stt_i, stm_t, stm_t_i, stt_t_i, t, t_s, tcm_time, tcm_idx, P_k_minus, P_k_plus, event_is_tcm, simparams)
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
maneuverSegments = simparams.maneuverSegments;
num_r_corrections = length(tcm_idx);

% Reshaping x so each column is a segment initial state and duration
x = reshape(x,m,n);

% Velocity mapping matrix
Mv = G';

% STM from beginning of trajectory to the state being targeted
% if simparams.target_final_maneuver
%     final_seg = simparams.maneuverSegments(end);
%     target_time = sum(x(7,1:final_seg - 1));
%     target_idx = find(target_time == t) ;
%     stmN0 = stm_t(:,:,target_idx);
%     stm_t = stm_t(:,:,1:target_idx);
%     t = t(1:target_idx);
%     t_s = t_s(1:target_idx);
% else
%     stmN0 = stm_t(:,:,end);
%     target_idx = length(t);
% end


% Store a gradient vector for each individual TCM (r and v)
% Structure for each tcm gradient
tcm_gradient_r = zeros(n*m, num_r_corrections);
tcm_gradient_v = zeros(n*m, 1);
% tcm_gradient_v_test = zeros(n*m, 1);

%% Loop through the corrections for each of the following sets of calcs

% Set up data structure for dPCkminusdxi: the partial derivative of the
% covariance matrix P (prior to correction C (the minus part)) with respect
% to each segment initial state
% For each correction there will be a 6x6x6 partial derivative of the
% covariance matrix P- with respect to each segment.
%%% to start - calculate for each correction and use each iteration, think
%%% there is no need to save all of them.

% Initialize partial derivative placeholders for each correction loop
dPCkminusdxi = zeros(6,6,6,length(event_is_tcm),n);
dPCkminusddti = zeros(6,6,length(event_is_tcm),n);
dTkdxi = zeros(3,6,6,length(event_is_tcm),n);
dTkddti = zeros(3,6,length(event_is_tcm),n);
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




%% New construct - looping through "events" instead of just TCMs
% - An event is a TCM or a node with a nominal maneuver (currently)





% event_times = zeros(num_events,1);
% event_times(event_is_tcm) = tcm_time;
% event_time_not_tcm = [];
% 
% for i = 1:sum(~event_is_tcm)
%     event_time_not_tcm(end+1) = sum(x(7,1:simparams.P_constrained_nodes(i) - 1));
% end
% 
% event_times(~event_is_tcm) = event_time_not_tcm';






[event_times, event_is_tcm] = define_events(x(:), t, tcm_time, simparams);


event_time_not_tcm = event_times(~event_is_tcm);
target_time = event_time_not_tcm(1);
target_idx = find(t==target_time);

num_events = length(event_times);


assert(length(event_times) == length(unique(event_times)),'Need some more logic to deal with a nominal maneuver and TCM simultaneous.')



for k = 1:num_events - 1

    
    % Get the time the kth event occurs
%     tCk = tcm_time(k);
    tCk = event_times(k);

    % Get the index tCk occurs
%     tCk_idx = tcm_idx(k);
    tCk_idx = find(t==tCk);

    % Test which segments have a correction
    tCk_seg = t_s(tCk_idx);
    tk_seg = t_s(tCk_idx);
    
    % Get the nominal state at the correction
%     x_tcm = x_t(tCk_idx,:)';

    % If the first iteration, the "last" TCM occured at t=0
    if k == 1
        tClast = 0;
        tClast_idx = 1;
        tlast = 0;
        tlast_idx = 1;
    end
    
    % Extract STMs for calculations    
    stmCkClast = dynCellCombine(t, t_s, tClast_idx, tCk_idx, simparams, stm_t_i);


    if event_times(k) == event_times(k+1)
        % A nominal maneuver and TCM are happening at the same time
        assert(0,'Need some more logic to deal with a nominal maneuver and TCM simultaneous.')
        

    elseif event_is_tcm(k)
        % Calculate the T matrix
        stmNCk = dynCellCombine(t, t_s, tCk_idx, target_idx, simparams, stm_t_i);
        Tk = [-inv( stmNCk(1:3,4:6) ) * stmNCk(1:3,1:3), -eye(3)];

        % Calculate the W and L matrices (used for TCMv calculation)
        Wk = [stmNCk(4:6,4:6), -stmNCk(4:6,1:3)];
        Lk = [Wk*Tk', zeros(3,3)];
    
        % Calculate N and IN matrices
        Nk = [zeros(3,6); Tk];
        INk = eye(6) + Nk;

    else
        % TODO: consider negating Tk so the wrong one doesn't get used, or something    
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

        
        [dstmCkClastdxi, dstmCkClastddti] = calc_dstmCkClast(x, x_i_f, t, t_s, stm_t, stm_i, stm_t_i, stt_t_i, tCk, tClast, i, simparams); % verified numerically



        %%%%%%% debugging/testing
%         final_target_time = event_times(end);
%         [dstmCkClastdxi, dstmCkClastddti] = calc_dstmCkClast(x, x_i_f, t, t_s, stm_t, stm_i, stm_t_i, stt_t_i, final_target_time, tCk, i, simparams); % verified numerically



        if event_is_tcm(k)
            % Assemble dTkdxi


            % dstmNC needs to change because the target is no longer at the
            % end of the trajectory, and that function assumes that it is.
            % Attempting to use calc_dstmCkClast instead:

            [dstmNCkdxi, dstmNCkddti] = calc_dstmCkClast(x, x_i_f, t, t_s, stm_t, stm_i, stm_t_i, stt_t_i, target_time, tCk, i, simparams);


%             [dstmNCkdxi, dstmNCkddti] = calc_dstmNC(x, x_i_f, t, t_s, stm_t, stm_i, stm_t_i, stt_t_i, tCk, tCk_idx, tCk_seg, i, simparams);
            dTkdxi(:,:,:,k,i) = T_partial(stmNCk, dstmNCkdxi);
            dTkddti(:,:,k,i) = T_partial(stmNCk, dstmNCkddti);
        end

        % Assemble dPCkminusdxi
        if k == 1 %%%% if the first correction
%             Tlast = zeros(3,6); dTlastdxi = zeros(3,6,6); dTlastddti = zeros(3,6); 
%             PClast_minus = simparams.P_initial; dPClast_minus_dxi = zeros(6,6); dPClast_minus_ddti = zeros(6,6);

            % Propagate dP to first event from P_initial
            dPCkminusdxi(:,:,:,k,i) = calc_dphiPphi(stmCkClast, dstmCkClastdxi, simparams.P_initial, zeros(6,6,6), zeros(6,6));
            dPCkminusddti(:,:,k,i) = calc_dphiPphi(stmCkClast, dstmCkClastddti, simparams.P_initial, zeros(6,6), zeros(6,6));

        else
            
            dTlastdxi = dTkdxi(:,:,:,k-1,i);
            dTlastddti = dTkddti(:,:,k-1,i);

            dPClast_minus_dxi = dPCkminusdxi(:,:,:,k-1,i);
            dPClast_minus_ddti = dPCkminusddti(:,:,k-1,i);

            % Structure P_k_minus has the dispersion P matrix prior to each EVENT 
            PClast_minus = P_k_minus(:,:,k-1);

            % Update differently if TCM vs. nominal maneuver event:
            % The event that matters is the last event (number k-1) because the
            % update to the current event time either involves a TCM or doesn't
            % pending what the previous event was.
            if event_is_tcm(k-1)
                dPCkminusdxi(:,:,:,k,i) = calc_dPckMinus(stmCkClast, dstmCkClastdxi, k, simparams, Tlast, dTlastdxi, PClast_minus, dPClast_minus_dxi);
                dPCkminusddti(:,:,k,i) = calc_dPckMinus(stmCkClast, dstmCkClastddti, k, simparams, Tlast, dTlastddti, PClast_minus, dPClast_minus_ddti);
    
            else   
                % Adding the impact of nominal maneuver execution error
                dPCkminusdxi(:,:,:,k,i) = calc_dphiPphi(stmCkClast, dstmCkClastdxi, PClast_minus, dPClast_minus_dxi, R_dv);
                dPCkminusddti(:,:,k,i) = calc_dphiPphi(stmCkClast, dstmCkClastddti, PClast_minus, dPClast_minus_ddti, R_dv);    
            end
                
            
        end       


        if event_is_tcm(k)
            % Assemble gradient
            % dxi
            tcm_gradient_r((i-1)*7+1:(i-1)*7+6,k) = calc_dSigk(Tk, dTkdxi(:,:,:,k,i), P_k_minus(:,:,k), dPCkminusdxi(:,:,:,k,i));
    
            % ddti
            tcm_gradient_r(i*7,k) = calc_dSigk(Tk, dTkddti(:,:,k,i), P_k_minus(:,:,k), dPCkminusddti(:,:,k,i));
        end

        % The final tcmV gradient
        if k == num_events - 1
            %% new method

            

            Pn = P_k_minus(:,:,end);
            dPndxi = calc_dPckMinus(stmNCk, dstmNCkdxi, k, simparams, Tk, dTkdxi(:,:,:,k,i), P_k_minus(:,:,k), dPCkminusdxi(:,:,:,k,i));
            dPnddti = calc_dPckMinus(stmNCk, dstmNCkddti, k, simparams, Tk, dTkddti(:,:,k,i), P_k_minus(:,:,k), dPCkminusddti(:,:,k,i));


            assert(event_is_tcm(k-1),'Error: need to incoroprate logic here for when the event prior to the final target is a nominal maneuver!');
      

            tcm_gradient_v((i-1)*7+1:(i-1)*7+6,1) = calc_dSigk(Mv, zeros(3,6,6), Pn, dPndxi);
            tcm_gradient_v(i*7,1) = calc_dSigk(Mv, zeros(3,6), Pn, dPnddti);
       



        end

        
    end % end of for loop for each segment (inside the loop for each TCM k)


    tClast = tCk; % Updating tClast for the next iteration
    tClast_idx = tCk_idx;
    tClast_seg = tCk_seg;
    


    % Assign new target time
    if event_is_tcm(k)
        % Updating Tlast, which is referenced the next iteration:
        % ( event_is_tcm(k-1) )
        Tlast = Tk;
    else
        if k < num_events
            target_time = event_time_not_tcm(find(event_time_not_tcm==target_time)+1);
            target_idx = find(t==target_time);     
        end
    end

end % end of loop for each TCM


%% Calculate output
tcm_gradient = sum(tcm_gradient_r')' + tcm_gradient_v;


end