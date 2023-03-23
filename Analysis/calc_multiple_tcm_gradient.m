function [tcm_gradient, tcm_gradient_r, tcm_gradient_v] = calc_multiple_tcm_gradient(x, x_t, x_i_f, stm_i, stt_i, stm_t, stm_t_i, stt_t_i, t, t_s, tcm_time, tcm_idx, tcm_dv_each, P_k_minus, P_k_plus, simparams)
%calc_tcm_gradient Computes and returns the analytical tcm gradient


%   Detailed explanation goes here

% Variable extraction
m = simparams.m;
n = simparams.n;
mu = simparams.mu;
P_initial = simparams.P_initial;
R = simparams.R;
G = [zeros(3,3); eye(3,3)];
maneuverSegments = simparams.maneuverSegments;
num_r_corrections = length(tcm_idx);

% Reshaping x so each column is a segment initial state and duration
x = reshape(x,m,n);

% Symplectic unit matrix
% J = [zeros(3,3), eye(3,3); -eye(3,3), zeros(3,3)];

% Velocity mapping matrix
Mv = [zeros(3,3), eye(3,3)];

% STM from beginning of trajectory to the state being targeted
if simparams.target_final_maneuver
    final_seg = simparams.maneuverSegments(end);
    target_time = sum(x(7,1:final_seg - 1));
    target_idx = find(target_time == t) ;
    stmN0 = stm_t(:,:,target_idx);
    stm_t = stm_t(:,:,1:target_idx);
    t = t(1:target_idx);
    t_s = t_s(1:target_idx);
else
    stmN0 = stm_t(:,:,end);
    target_idx = length(t);
end


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
dPCkminusdxi = zeros(6,6,6,num_r_corrections,n);
dPCkminusddti = zeros(6,6,num_r_corrections,n);
dTkdxi = zeros(3,6,6,num_r_corrections,n);
dTkddti = zeros(3,6,num_r_corrections,n);
% dLkdxi = zeros(3,6,6,num_r_corrections,n);
% dLkddti = zeros(3,6,num_r_corrections,n);
% dstmNCkdxi = zeros(6,6,6,num_r_corrections,n);
% dstmNCkddti = zeros(6,6,num_r_corrections,n);
% dPndxi = zeros(6,6,6,n);
% dPddti = zeros(6,6,n);

for k = 1:num_r_corrections

    
    

    % Get the time the kth correction occurs
    tCk = tcm_time(k);

    % Get the index tCk occurs
    tCk_idx = find(t==tCk);

    % Test which segments have a correction
    tCk_seg = t_s(tCk_idx);
    
    % Get the nominal state at the correction
    x_tcm = x_t(tCk_idx,:)';

    % If the first iteration, the "last" TCM occured at t=0
    if k == 1
        tClast = 0;
        tClast_idx = 1;
    end
    
    % Extract STMs for calculations    
    stmCkClast = dynCellCombine(t, t_s, tClast_idx, tCk_idx, simparams, stm_t_i);
    stmNCk = dynCellCombine(t, t_s, tCk_idx, target_idx, simparams, stm_t_i);

    % Calculate the T matrix
    Tk = [-inv( stmNCk(1:3,4:6) ) * stmNCk(1:3,1:3), -eye(3)];
    
    % Calculate the W and L matrices (used for TCMv calculation)
    Wk = [stmNCk(4:6,4:6), -stmNCk(4:6,1:3)];
    Lk = [Wk*Tk', zeros(3,3)];

    % Calculate N and IN matrices
    Nk = [zeros(3,6); Tk];
    INk = eye(6) + Nk;
    
    
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

        [dstmNCkdxi, dstmNCkddti] = calc_dstmNC(x, x_i_f, t, t_s, stm_t, stm_i, stm_t_i, stt_t_i, tCk, tCk_idx, tCk_seg, stmN0, i, simparams);


        % Assemble dTkdxi
        dTkdxi(:,:,:,k,i) = T_partial(stmNCk, dstmNCkdxi);
        dTkddti(:,:,k,i) = T_partial(stmNCk, dstmNCkddti);

        % Assemble dPCkminusdxi
        if k == 1 %%%% if the first correction
            Tlast = zeros(3,6); dTlastdxi = zeros(3,6,6); dTlastddti = zeros(3,6); 
            PClast_minus = zeros(6,6); dPClast_minus_dxi = zeros(6,6); dPClast_minus_ddti = zeros(6,6);
        else
            
            dTlastdxi = dTkdxi(:,:,:,k-1,i);
            dTlastddti = dTkddti(:,:,k-1,i);

            dPClast_minus_dxi = dPCkminusdxi(:,:,:,k-1,i);
            dPClast_minus_ddti = dPCkminusddti(:,:,k-1,i);

            % Structure P_i_minus has the dispersion P matrix prior to each correction
            PClast_minus = P_k_minus(:,:,k-1);
            
        end

        
        dPCkminusdxi(:,:,:,k,i) = calc_dPckMinus(stmCkClast, dstmCkClastdxi, k, simparams, Tlast, dTlastdxi, PClast_minus, dPClast_minus_dxi);
        dPCkminusddti(:,:,k,i) = calc_dPckMinus(stmCkClast, dstmCkClastddti, k, simparams, Tlast, dTlastddti, PClast_minus, dPClast_minus_ddti);

        % Assemble gradient
        % dxi
        tcm_gradient_r((i-1)*7+1:(i-1)*7+6,k) = calc_dSigk(Tk, dTkdxi(:,:,:,k,i), P_k_minus(:,:,k), dPCkminusdxi(:,:,:,k,i));

        % ddti
        tcm_gradient_r(i*7,k) = calc_dSigk(Tk, dTkddti(:,:,k,i), P_k_minus(:,:,k), dPCkminusddti(:,:,k,i));




        if k == num_r_corrections
            %% new method

            Pn = P_k_minus(:,:,k+1);
            dPndxi = calc_dPckMinus(stmNCk, dstmNCkdxi, k, simparams, Tk, dTkdxi(:,:,:,k,i), P_k_minus(:,:,k), dPCkminusdxi(:,:,:,k,i));
            dPnddti = calc_dPckMinus(stmNCk, dstmNCkddti, k, simparams, Tk, dTkddti(:,:,k,i), P_k_minus(:,:,k), dPCkminusddti(:,:,k,i));
      

            tcm_gradient_v((i-1)*7+1:(i-1)*7+6,1) = calc_dSigk(Mv, zeros(3,6,6), Pn, dPndxi);
            tcm_gradient_v(i*7,1) = calc_dSigk(Mv, zeros(3,6), Pn, dPnddti);
       


        end

        
    end % end of for loop for each segment (inside the loop for each TCM k)


    tClast = tCk; % Updating tClast for the next iteration
    tClast_idx = tCk_idx;
    tClast_seg = tCk_seg;
    % TODO: ALSO ASSIGN LASTS FOR Tlast, dTlastdxi, PClast_minus, AND dPClast_minus_dxi
    Tlast = Tk;

end % end of loop for each TCM

%% Calculate tcm_gradient_v

% TCM_V calculations

% clear('dPCkminusdxi','dPCkminusddti','dTkdxi','dTkddti')
% 
% for i = 1:n
% 
% 
% 
%     Tlast = zeros(3,6); dTlastdxi = zeros(3,6,6); dTlastddti = zeros(3,6); 
% %     PClast_minus = zeros(6,6); 
%     PClast_minus = simparams.P_initial; 
%     dPCkminusdxi = zeros(6,6,6); dPCkminusddti = zeros(6,6);
%     tClast = 0;
%     tClast_idx = 1;
% 
%     for k = 1:num_r_corrections
% 
% 
%         if i == 19
%             ppp=1; % debug spot
%         end
% 
%         % Get the time the kth correction occurs
%         tCk = tcm_time(k);    
%         % Get the index tCk occurs
%         tCk_idx = find(t==tCk);    
%         % Test which segments have a correction
%         tCk_seg = t_s(tCk_idx);
%         
%         % Extract STMs for calculations    
%         stmCkClast = dynCellCombine(t, t_s, tClast_idx, tCk_idx, simparams, stm_t_i);
%         stmNCk = dynCellCombine(t, t_s, tCk_idx, target_idx, simparams, stm_t_i);
%     
%         % Calculate the T matrix
%         Tk = [-inv( stmNCk(1:3,4:6) ) * stmNCk(1:3,1:3), -eye(3)];
% 
%         % STM derivatives
%         [dstmCkClastdxi, dstmCkClastddti] = calc_dstmCkClast(x, x_i_f, t, t_s, stm_t, stm_i, stm_t_i, stt_t_i, tCk, tClast, i, simparams); % verified numerically
%         [dstmNCkdxi, dstmNCkddti] = calc_dstmNC(x, x_i_f, t, t_s, stm_t, stm_i, stm_t_i, stt_t_i, tCk, tCk_idx, tCk_seg, stmN0, i, simparams);
% 
% %         % Propagate P 
% %         if i == 3 && k == 1
% %             pppp=1;
% %         end
%         dPCkminusdxi = calc_dPckMinus(stmCkClast, dstmCkClastdxi, k, simparams, Tlast, dTlastdxi, PClast_minus, dPCkminusdxi);
%         dPCkminusddti = calc_dPckMinus(stmCkClast, dstmCkClastddti, k, simparams, Tlast, dTlastddti, PClast_minus, dPCkminusddti);
% 
% 
% %         % DEBUGGING%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %         % Print out the finite difference equivalent for comparison
% %         i
% %         k
% %         
% % 
% %         if k == 1
% %             dPCkminusdxi_fd = dPc1_minusdx_fd(:,:,(i-1)*7+1:(i-1)*7+6);
% %             
% % 
% %             dPCkminusddti_fd = dPc1_minusdx_fd(:,:,(i-1)*7+7);
% %             
% %         elseif k == 2
% %             dPCkminusdxi_fd = dPc2_minusdx_fd(:,:,(i-1)*7+1:(i-1)*7+6);
% % 
% %             dPCkminusddti_fd = dPc2_minusdx_fd(:,:,(i-1)*7+7);
% %             if i == 2
% %                 ppp=1; % debug spot
% %             end
% %         elseif k == 3
% %             dPCkminusdxi_fd = dPc3_minusdx_fd(:,:,(i-1)*7+1:(i-1)*7+6);
% % 
% %             dPCkminusddti_fd = dPc3_minusdx_fd(:,:,(i-1)*7+7);
% %         elseif k == 4
% %             dPCkminusdxi_fd = dPc4_minusdx_fd(:,:,(i-1)*7+1:(i-1)*7+6);
% % 
% %             dPCkminusddti_fd = dPc4_minusdx_fd(:,:,(i-1)*7+7);
% %         end
% % 
% %         (dPCkminusdxi_fd - dPCkminusdxi)./dPCkminusdxi_fd
% %         (dPCkminusddti_fd - dPCkminusddti)./dPCkminusddti_fd
% %         
% % 
% %         % DEBUGGING END%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%         %%% The following for the next iteration %%%
%         % Assemble dTkdxi
%         dTlastdxi = T_partial(stmNCk, dstmNCkdxi);
%         dTlastddti = T_partial(stmNCk, dstmNCkddti);
% 
%         tClast = tCk; % Updating tClast for the next iteration
%         tClast_idx = tCk_idx;
%         tClast_seg = tCk_seg;
%         Tlast = Tk;
%         PClast_minus = P_k_minus(:,:,k);
% 
%     end
% 
% 
%     Pn = P_k_minus(:,:,k+1);
% 
%     dPndxi = calc_dPckMinus(stmNCk, dstmNCkdxi, k, simparams, Tk, dTlastdxi, P_k_minus(:,:,k), dPCkminusdxi);
%     dPnddti = calc_dPckMinus(stmNCk, dstmNCkddti, k, simparams, Tk, dTlastddti, P_k_minus(:,:,k), dPCkminusddti);
% 
% 
% 
% % 
% %     % DEBUGGING%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     % Print out the finite difference equivalent for comparison
% %     i
% %     k
% %         dPCnminusdxi_fd = dPcn_minusdx_fd(:,:,(i-1)*7+1:(i-1)*7+6);
% %         (dPCnminusdxi_fd - dPndxi)./dPCnminusdxi_fd
% % 
% %         dPCnminusddti_fd = dPc1_minusdx_fd(:,:,(i-1)*7+7);
% %         (dPCnminusddti_fd - dPnddti)./dPCnminusddti_fd
% %  
% %     
% % 
% %     % DEBUGGING END%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % 
% 
% 
% 
% 
% %     tcm_gradient_v((i-1)*7+1:(i-1)*7+6,1) = calc_dSigk(Mv, zeros(3,6,6), Pn, dPndxi);
% %     tcm_gradient_v(i*7,1) = calc_dSigk(Mv, zeros(3,6), Pn, dPnddti);
%     
% 
% end



tcm_gradient = sum(tcm_gradient_r')' + tcm_gradient_v;


end