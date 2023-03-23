function [tcm_dv_total] = calc_tcmdv(x, t, stm_t, tcm_idx, simparams)
%calc_covariance_tcmdv Calculates the final dispersion covariance and the
%total delta V of the position corrections occuring along the nominal
%trajectory at each entry in the array tcm_time, as well as calculating the
%target dispersion covariance (P_target)



tcm_time = t(tcm_idx)';

P_initial = simparams.P_initial;

m = simparams.m; % the number of elements per segment
n = simparams.n; % the number of segments

% Reshaping x so each column is a segment initial state and duration
x = reshape(x,m,n);


R = zeros(6,6);
R(4:6,4:6) = simparams.R; % TCM execution error covariance

% STM from beginning of trajectory to the state being targeted

if simparams.target_final_maneuver
    final_seg = simparams.maneuverSegments(end);
    target_time = sum(x(7,1:final_seg - 1));
    target_idx = find(target_time == t) ;
    stmN0 = stm_t(:,:,target_idx);

else
    stmN0 = stm_t(:,:,end);
end




% Symplectic unit matrix
J = [zeros(3,3), eye(3,3); -eye(3,3), zeros(3,3)];


%%%%%% testing functionality if add a second maneuver
%%%%%%%  REMOVE THIS PART AFTERWARDS

% tcm_time = [tcm_time, t(318)];
% tcm_time = [t(318)];

num_tcm = length(tcm_time);
%%%%%% REMOVE TO HERE

P_i = P_initial;

P_t = zeros(6,6,length(t));
% Empty structure to store the tcm_dv RSS magnitudes (not 3 sigma). It is
% one longer than the tcm_time because the velocity correction occurs at
% the end of the trajectory and isn't included in tcm_time.
tcm_dv = zeros(1,length(tcm_time)+1);


%% Logic if no TCM (zero length)
if length(tcm_time) == 0

    P_target = stmN0 * P_i * stmN0';
    tcm_dv_total = 0;

else % Otherwise, if there are TCMs, do:

    %% Setup for logic choices and STM portions
    % Test which segment has the correction
    
    tcm_idx_logical = logical(sum(t'==tcm_time', 1));
    tcm_t_idx = find(tcm_idx_logical);

    if sum(tcm_idx_logical)~=num_tcm
        assert(0,'Error: you likely input duplicate TCM times in your tcm_time array');
    end
    
    % corrSeg = t_s(tcm_idx_logical);
    
    % % Get the nominal state at the correction (DON'T THINK IS NEEDED, BUT
    % SAVING JUST IN CASE)
    % x_tcm = x_t(tcm_idx_logical,:)';
    
    % Extract STMs for calculations
    % A tensor of STMs from the beginning of the trajectory to each correction
    stmC0 = stm_t(:,:,tcm_idx_logical);
    



        
    
    
    %% Compute the covariance only at the correction(s) and the final
    for i = 1:num_tcm
        % Propagate dispersion covariance from previous event to i
        if i == 1
            stmCiClast = stmC0(:,:,i);
        else
            stmCi0 = stmC0(:,:,i);
    %         stmClast0 = stmC0(:,:,i-1);
    %         stm0Clast = -J * stmClast0' * J;
            stmCiClast = stmCi0 * stm0C;
        end
    
        %%%%%%%%%%%%% DEBUGGING %%%%%%%%%%%%%%%%%%
        % Add a test  - separately propagate stmCiClast to verify %%%%%%%%%%%%
    %     if i == 1
    %         xi = x(1:6)';
    %     else
    %         xi = x_t(tcm_t_idx(i-1),:)';
    %     end
    %     xi = x(1:6)';
    %     dt = tcm_time(i);
    %     
    % %     [x2, stm21] = stateStmProp(xi,dt,simparams);
    % 
    % 
    %     [x2,stm1] = stateStmProp(xi,x(7,1),simparams);
    %     [x3,stm2] = stateStmProp(x(1:6,2),x(7,2),simparams);
    % 
    %     [xc2, stmc22] = stateStmProp(x(1:6,2),tcm_time(2) - x(7,1),simparams)
    % 
    % 
    %     % also try three (?) propagations, to get from tcm1 to tcm2
    %     % get state at tcm1
    %     [xtcm1, stmc11] = stateStmProp(xi,tcm_time(1),simparams);
    %     % propagate to end of seg 1 to get the stm
    %     [x2,stm2c1] = stateStmProp(xtcm1,x(7,1)-tcm_time(1),simparams);
    %     stmc22 * stm2c1
    % 
    %     % get stmNC
    %     % Pro
    %     [x3,stmNc2] = stateStmProp(xc2,x(7,1)+x(7,2)-tcm_time(2),simparams);
    % 
    %     % stmN0
    %     stmNc2 * stmc22 * stm1
    % 
    % 
    %     stmc22 * stm2
    % 
    % %     stmf
    %     stmC0(:,:,i)
    % 
    %     % Separately do the whole thing
    % 
    %     Pc1_minus = stmc11 * P_initial * stmc11'
    %     stmNc1 = stm2 * stm2c1;
    % 
    %     T = [-inv( stmNc1(1:3,4:6) ) * stmNc1(1:3,1:3), -eye(3)];
    %     N = [zeros(3,6); T];
    %     IN = eye(6) + N;
    % 
    %     Pc1_plus = IN*Pc1_minus * IN' + R
    %     stmc2c1 = stmc22*stm2c1
    %     Pc2_minus = stmc2c1 * Pc1_plus * stmc2c1'
    % 
    %     T = [-inv( stmNc2(1:3,4:6) ) * stmNc2(1:3,1:3), -eye(3)];
    %     N = [zeros(3,6); T];
    %     IN = eye(6) + N;
    % 
    % 
    %     Pc2_plus = IN * Pc2_minus * IN' + R
    %     Pf = stmNc2 * Pc2_plus * stmNc2'
    
    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        
        P_i_minus = stmCiClast * P_i * stmCiClast';
    
        % Perform the TCM update
        stm0C = -J * stmC0(:,:,i)' * J;
        stmNC = stmN0 * stm0C;
    
        T = [-inv( stmNC(1:3,4:6) ) * stmNC(1:3,1:3), -eye(3)];
        N = [zeros(3,6); T];
        IN = eye(6) + N;
    
    
    %     if i == 2
    %         R = zeros(6,6);
    %     end
    
        P_i = IN * P_i_minus * IN' + R;

        P_tcm = T * P_i_minus * T' + simparams.R;
        tcm_dv(i) = sqrt(trace(P_tcm));
        
    end
    
    % Propagate the rest of the way to the target
    % stmClast0 = stmC0(:,:,end);
    % stm0Clast = -J * stmClast0' * J;
    % stmNClast = stmN0 * stm0Clast;
    
    P_target = stmNC * P_i * stmNC';

    tcm_dv(end) = sqrt(trace(P_target(4:6,4:6)));

    tcm_dv_total = sum(tcm_dv);

end


end