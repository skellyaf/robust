function [Q_k_km1, dQ_k_km1_dxi, dQ_k_km1_ddti] = calc_Q_events(traj, x, tcm_time, simparams)
%calc_Q_events returns the Qbar from the previous event to the current
%event, before any covariance update occurs at the event. The covariance
%update then resets Qbar and it starts again at zero after each event.

x = reshape(x,simparams.m,simparams.n);
%% Define events

% Get event info
[event_times, event_indicator] = define_events_v2(x(:), traj.t, tcm_time, simparams);

event_idx_logical = logical(sum(traj.t'==event_times', 1));    
event_idxs = find(event_idx_logical);

Q_k_km1 = zeros(6,6,length(event_indicator));

if nargout > 1
    dQ_k_km1_dxi = zeros(6,6,6,length(event_indicator),simparams.n);

    dQ_k_km1_ddti = zeros(6,6,length(event_indicator),simparams.n);
end


if simparams.start_P_growth_node > 1
    t_km1 = sum(x(7,1:simparams.start_P_growth_node-1));
    idx_km1 = find(traj.t == t_km1);
elseif simparams.start_P_growth_node == 1
    idx_km1 = 1;
    t_km1 = 0;
end

for k = 1:length(event_indicator) 
    % If it a dispersion covariance modifying event
    if event_indicator(k) > 0 || k == length(event_indicator)

        % What is the current event index
        idx_k = event_idxs(k);
        t_k = event_times(k);

        % Calculate accumulated Q from the previous event to the current event
        Q_k_km1(:,:,k) = Qcombine(traj, idx_km1, idx_k, simparams);

        % Calculate partial wrt each segment initial state, loop through
        % segments
        if nargout > 1
            for i = 1:simparams.n
                if i == 1
                    t0_i = 0;
                    idx_i = 1;
                else
                    t0_i = sum(x(7,1:i-1));
                    idx_i = find(traj.t == t0_i);
                end

                tf_i = sum(x(7,1:i));
                idx_fi = find(traj.t == tf_i);

                if isempty(idx_fi)
                    [~, idx_fi] = min(abs(traj.t-tf_i));
                end
                    




                %%%%% Qbar state sensitivity
                if idx_k < idx_i || idx_fi < idx_km1 % attempt with indices to fix a bug
%                 if t_k < t0_i || tf_i < t_km1
                    % Already populated with zeros, do nothing.
%                     dQ_k_km1_dxi(:,:,:, k, i) = zeros(6,6,6);

                elseif idx_i <= idx_km1 && idx_fi <= idx_k % CASE 1 %%%%% NUMERICALLY VERIFIED
%                 elseif t0_i <= t_km1 && tf_i <= t_k % CASE 1 %%%%% NUMERICALLY VERIFIED
                    
                    stm_k_fi = dynCellCombine(traj.t, traj.t_s, idx_fi, idx_k, simparams, traj.stm_t_i);

                    stm_fi_km1 = dynCellCombine(traj.t, traj.t_s, idx_km1, idx_fi, simparams, traj.stm_t_i);
                    dstm_fi_km1_dx0i = calc_dstmCkClast(x, traj.x_i_f, traj.t, traj.t_s, traj.stm_t, traj.stm_i, traj.stm_t_i, traj.stt_t_i, tf_i, t_km1, i, simparams);
                    
                    dQ_fi_0i_dx0i = traj.dQ_t_i{i}(:,:,:,end);
                    
                    km1_i_idx = find_cell_t_i_idx_2(traj.t, traj.t_s, idx_km1, traj.stm_t_i, i);
                    Q_km1_0i = traj.Q_t_i{i}(:,:,km1_i_idx);
                    dQ_km1_0i_dx0i = traj.dQ_t_i{i}(:,:,:,km1_i_idx);


                    dQ_fi_km1_dx0i = dQ_fi_0i_dx0i ...
                                   - tmult(dstm_fi_km1_dx0i,tmult(Q_km1_0i,stm_fi_km1,[0 1])) ...
                                   - tmult(stm_fi_km1,tmult(dQ_km1_0i_dx0i,stm_fi_km1,[0 1])) ...
                                   - tmult(stm_fi_km1,tmult(Q_km1_0i,dstm_fi_km1_dx0i,[0 1]));

                    

                    
                    dQ_k_km1_dxi(:,:,:,k,i) = tmult(stm_k_fi, tmult(dQ_fi_km1_dx0i, stm_k_fi, [0,1]));

                    % dt

% % %                     if tf_i == t_k
% % %                         x_k = traj.x_t(idx_k,:)';
% % %                         dQ_k_km1_ddti(:,:,k,i) = Qdot(Q_k_km1(:,:,k), x_k, simparams.mu, simparams.Qt);
% % % 
% % % 
% % %                     end

%                     [Q_fi_km1, dQ_fi_km1] = Qcombine(traj, idx_km1, idx_fi, simparams);
%                     x_fi = traj.x_t(idx_fi,:)';
                    
%                     dQ_fi_km1_ddti = Qdot(Q_fi_km1, x_fi, simparams.mu, simparams.Qt);
%                     dQ_k_km1_ddti(:,:,k,i) = stm_k_fi * dQ_fi_km1_ddti * stm_k_fi'; 




                elseif idx_i <= idx_km1 && idx_k < idx_fi % CASE 2 %%%%% NUMERICALLY VERIFIED
%                 elseif t0_i <= t_km1 && t_k < tf_i % CASE 2 %%%%% NUMERICALLY VERIFIED
                    

                    k_i_idx = find_cell_t_i_idx_2(traj.t,traj.t_s,idx_k, traj.stm_t_i, i);
                    km1_i_idx = find_cell_t_i_idx_2(traj.t,traj.t_s,idx_km1, traj.stm_t_i, i);

                    Q_km1_0i = traj.Q_t_i{i}(:,:,km1_i_idx);
                    dQ_k_0i_dx0i = traj.dQ_t_i{i}(:,:,:,k_i_idx);
                    dQ_km1_0i_dx0i = traj.dQ_t_i{i}(:,:,:,km1_i_idx);

                    dstm_k_km1_dx0i = calc_dstmCkClast(x, traj.x_i_f, traj.t, traj.t_s, traj.stm_t, traj.stm_i, traj.stm_t_i, traj.stt_t_i, t_k, t_km1, i, simparams);

                    stm_k_km1 = dynCellCombine(traj.t, traj.t_s, idx_km1, idx_k, simparams, traj.stm_t_i);

                    dQ_k_km1_dxi(:,:,:,k,i) = dQ_k_0i_dx0i ...
                                            - tmult(dstm_k_km1_dx0i,tmult(Q_km1_0i,stm_k_km1,[0 1])) ...
                                            - tmult(stm_k_km1,tmult(dQ_km1_0i_dx0i,stm_k_km1,[0 1])) ...
                                            - tmult(stm_k_km1,tmult(Q_km1_0i,dstm_k_km1_dx0i,[0 1]));

                    % dt is zero unless t_k = tf_i, but that situation is
                    % covered by CASE 1
                    

                elseif idx_km1 < idx_i && idx_fi <= idx_k % CASE 3 %%%%% NUMERICALLY VERIFIED
%                 elseif t_km1 < t0_i && tf_i <= t_k % CASE 3 %%%%% NUMERICALLY VERIFIED

                    % dQk_tfi_t0i_dx0i
                    dQ_fi_0i_dx0i = traj.dQ_t_i{i}(:,:,:,end);
                    stm_k_fi = dynCellCombine(traj.t, traj.t_s, idx_fi, idx_k, simparams, traj.stm_t_i);
                    dQk_fi_0i_dx0i = tmult(stm_k_fi,tmult(dQ_fi_0i_dx0i,stm_k_fi,[0,1]));

                    % dQk_t0i_tkm1_dx0i
                    stm_fi_0i = traj.stm_t_i{i}(:,:,end);
                    stt_fi_0i = traj.stt_t_i{i}(:,:,:,end);
                    Q_0i_km1 = Qcombine(traj, idx_km1, idx_i, simparams);
                    dQfi_0i_km1_dx0i = tmult(stt_fi_0i, tmult(Q_0i_km1, stm_fi_0i, [0,1])) ...
                                     + tmult(stm_fi_0i, tmult(Q_0i_km1, stt_fi_0i, [0,1]));

                    dQk_0i_km1_dx0i = tmult(stm_k_fi, tmult(dQfi_0i_km1_dx0i, stm_k_fi, [0,1]));


                    dQ_k_km1_dxi(:,:,:,k,i) = dQk_fi_0i_dx0i + dQk_0i_km1_dx0i;



                    % dt
%                     x_k = traj.x_t(idx_k,:)';
%                     dQ_k_km1_ddti(:,:,k,i) = Qdot(Q_k_km1(:,:,k), x_k, simparams.mu, simparams.Qt);



                    %%%% test
% % % % %                     Q_fi_0i = traj.Q_t_i{i}(:,:,end);
% % % % %                     x_fi = traj.x_t(idx_fi,:)';
% % % % % 
% % % % %                     dQ_fi_0i_ddti = Qdot(Q_fi_0i, x_fi, simparams.mu, simparams.Qt);
% % % % % 
% % % % % %                     stm_k_fi
% % % % % 
% % % % %                     Q_k_fi = Qcombine(traj, idx_fi, idx_k, simparams);
% % % % %                     dQ_k_fi_ddti = Qdot(Q_k_fi, x_k, simparams.mu, simparams.Qt);
% % % % % 
% % % % % %                     Q_0i_km1 
% % % % %                     x0_i = x(1:6,i);
% % % % %                     dQ_0i_km1_ddti = Qdot(Q_0i_km1, x0_i, simparams.mu, simparams.Qt);
% % % % %                     stm_k_0i = dynCellCombine(traj.t, traj.t_s, idx_i, idx_k, simparams, traj.stm_t_i);
% % % % % 
% % % % % 
% % % % % 
% % % % % %                     tes = stm_k_0i * dQ_0i_km1_ddti * stm_k_0i' ...
% % % % % %                         + stm_k_fi * dQ_fi_0i_ddti * stm_k_fi'
% % % % % %                         - dQ_k_fi_ddti
% % % % % 
% % % % %                     stmdot_
% % % % %                     tes2 = stm_k_fi * dQ_fi_0i_ddti * stm_k_fi' ...
% % % % %                         - dQ_k_fi_ddti 
% % % % %                    
% % % % % 
% % % % % 
% % % % %                     %%%% test end
% % % % % 
% % % % %                     Q_fi_0i = traj.Q_t_i{i}(:,:,end);
% % % % %                     x_fi = traj.x_t(idx_fi,:)';
% % % % % 
% % % % %                     dQ_fi_0i_ddti = Qdot(Q_fi_0i, x_fi, simparams.mu, simparams.Qt);
% % % % %                     dQk_fi_0i_ddti = tmult(stm_k_fi, tmult(dQ_fi_0i_ddti ,stm_k_fi ,[0,1])); %%%%%% TODO -REGULAR MATRIX MULTIPLICATION IS FINE
% % % % % 
% % % % % 
% % % % %                     stmdot_fi = stmDot(x_fi, stm_fi_0i, simparams);
% % % % %                     dQk_0i_km1_ddti = stm_k_fi * stmdot_fi * Q_0i_km1 * stm_fi_0i' * stm_k_fi' ...
% % % % %                                     + stm_k_fi * stm_fi_0i * Q_0i_km1 * stmdot_fi' * stm_k_fi';
% % % % % 
% % % % % 
% % % % %                     dQ_k_km1_ddti(:,:,k,i) = dQk_fi_0i_ddti + dQk_0i_km1_ddti;

                    




                elseif idx_km1 < idx_i && idx_k < idx_fi % CASE 4 %%%%% NUMERICALLY VERIFIED
%                 elseif t_km1 < t0_i && t_k < tf_i % CASE 4 %%%%% NUMERICALLY VERIFIED

                    k_i_idx = find_cell_t_i_idx_2(traj.t,traj.t_s,idx_k, traj.stm_t_i, i);

                    dQ_k_0i_dx0i = traj.dQ_t_i{i}(:,:,:,k_i_idx);

                    Q_0i_km1 = Qcombine(traj, idx_km1, idx_i, simparams);
                    stm_k_0i = traj.stm_t_i{i}(:,:,k_i_idx);
                    stt_k_0i = traj.stt_t_i{i}(:,:,:,k_i_idx);
                    dQk_0i_km1_dx0i = tmult(stt_k_0i, tmult(Q_0i_km1, stm_k_0i, [0,1])) ...
                                    + tmult(stm_k_0i, tmult(Q_0i_km1, stt_k_0i, [0,1]));

                    dQ_k_km1_dxi(:,:,:,k,i) = dQ_k_0i_dx0i + dQk_0i_km1_dx0i;

                    % dt is nonzero when t_k equals tf_i, which is taken
                    % care of by the case above!

                    
                    
                    


                else
                    assert(0,'Something went wrong!')
                end


                %%%%% Sensitivity to segment duration
                if nargout > 1
                    if idx_km1 < idx_fi && idx_fi <= idx_k
%                     if t_km1 < tf_i && tf_i <= t_k
                        Q_fi_km1 = Qcombine(traj, idx_km1, idx_fi, simparams);
                        x_fi = traj.x_t(idx_fi,:)';
                        stm_k_fi = dynCellCombine(traj.t, traj.t_s, idx_fi, idx_k, simparams, traj.stm_t_i);
                        dQ_fi_km1_ddti = Qdot(Q_fi_km1, x_fi, simparams.mu, simparams.Qt);
                        dQ_k_km1_ddti(:,:,k,i) = stm_k_fi * dQ_fi_km1_ddti * stm_k_fi';




        
                    end
                end




            end
        end

        % After calcs are done, update idx_prev_event for next round
        idx_km1 = event_idxs(k);
        t_km1 = event_times(k);
    end
end












end