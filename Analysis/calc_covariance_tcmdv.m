function [P_target, tcm_dv_total, tcm_dv, P_i_minus, P_i_plus] = calc_covariance_tcmdv(t, stm_t, tcm_time, vel_disp_flag, P_i, simparams)
%calc_covariance_tcmdv Calculates the final dispersion covariance and the
%total delta V of the position corrections occuring along the nominal
%trajectory at each entry in the array tcm_time, as well as calculating the
%target dispersion covariance (P_target)




% % OPTIONAL INPUT: P_i
% if length(varargin) < 1
%     P_i = simparams.P_initial;
% else
%     P_i = varargin{1};
% end


% m = simparams.m; % the number of elements per segment
% n = simparams.n; % the number of segments
% 
% % Reshaping x so each column is a segment initial state and duration
% x = reshape(x,m,n);

R = zeros(6,6);
R(4:6,4:6) = simparams.R; % TCM execution error covariance

% STM from beginning of trajectory to the state being targeted

% if simparams.target_final_maneuver
%     final_seg = simparams.maneuverSegments(end);
%     target_time = sum(x(7,1:final_seg - 1));
%     target_idx = find(target_time == t) ;
%     stmN0 = stm_t(:,:,target_idx);
% 
% else
%     stmN0 = stm_t(:,:,end);
% end
% 
% if length(t) ~= target_idx
%     ppp=1;
% end
stmN0 = stm_t(:,:,end);

% Empty structure to store the tcm_dv RSS magnitudes (not 3 sigma). It is
% one longer than the tcm_time because the velocity correction occurs at
% the end of the trajectory and isn't included in tcm_time.
num_tcm = length(tcm_time);
tcm_dv = zeros(1,length(tcm_time)+1);



%% Logic if no TCM (zero length)
% if length(tcm_time) == 0
if isempty(tcm_time)

%     METHOD 2 PART:
%     stmN0 = dynCellCombine(t, t_s, start_idx, target_idx, simparams, stm_t_i);

    P_target = stmN0 * P_i * stmN0';
    P_i_minus = P_target;
    P_i_plus = nan;
    tcm_dv_total = 0;

else % Otherwise, if there are TCMs, do:

    %% Setup for logic choices and STM portions
    % Test which segment has the correction
    
    tcm_idx_logical = logical(sum(t'==tcm_time', 1));
%     C_idxs = find(tcm_idx_logical);

    Clast_idx = 1;

    if sum(tcm_idx_logical)~=num_tcm
        assert(0,'Error: you likely input duplicate TCM times in your tcm_time array');
    end
    
    % corrSeg = t_s(tcm_idx_logical);
    
    % % Get the nominal state at the correction (DON'T THINK IS NEEDED, BUT
    % SAVING JUST IN CASE)
    % x_tcm = x_t(tcm_idx_logical,:)';
    
    % Extract STMs for calculations
    % A tensor of STMs from the beginning of the trajectory to each correction

    % METHOD 1 PART
    stmC0 = stm_t(:,:,tcm_idx_logical);

    
    


    % Empty structure for storing the before correction and after
    % correction dispersion covariance
    P_i_minus = zeros(6,6,num_tcm+1); % plus one because the last is going to be P at the target
    P_i_plus = zeros(6,6,num_tcm); % no plus one, it is the dispersion covariance after each tcm

        
    
    
    %% Compute the covariance only at the correction(s) and the final
    for i = 1:num_tcm
        % Propagate dispersion covariance from previous event to i

        % METHOD 1
        if i == 1
            stmCiClast = stmC0(:,:,i);
        else
            stmCi0 = stmC0(:,:,i);
            stmCiClast = stmCi0 * stm0C;
        end 

        % METHOD 2
%         C_idx = C_idxs(i);
%         stmCiClast = dynCellCombine_r2(t, t_s, Clast_idx, C_idx, simparams, stm_t_i);
        
        P_i_minus(:,:,i) = stmCiClast * P_i * stmCiClast';
    
        % Perform the TCM update


        % METHOD 1
        stm0C = invert_stm(stmC0(:,:,i), simparams);
        stmNC = stmN0 * stm0C;

        % METHOD 2
%         stmNC = dynCellCombine_r2(t, t_s, C_idx, target_idx, simparams, stm_t_i);


    
        T = [-inv( stmNC(1:3,4:6) ) * stmNC(1:3,1:3), -eye(3)];
        N = [zeros(3,6); T];
        IN = eye(6) + N;
   
        P_i = IN * P_i_minus(:,:,i) * IN' + R;
        P_i_plus(:,:,i) = P_i;

        P_tcm = T * P_i_minus(:,:,i) * T' + simparams.R; % with error in tcm_dv magnitude calc 
%         P_tcm = T * P_i_minus(:,:,i) * T'; % without error in tcm_dv magnitude calc 
        tcm_dv(i) = sqrt(trace(P_tcm));

        % METHOD 2 PART
%         Clast_idx = C_idx;
        
    end
    
    % Propagate the rest of the way to the target
    
    P_target = stmNC * P_i * stmNC';
    P_i_minus(:,:,end) = P_target;

    if vel_disp_flag
        tcm_dv(end) = sqrt(trace(P_target(4:6,4:6)));
    end

    tcm_dv_total = sum(tcm_dv);

end


end