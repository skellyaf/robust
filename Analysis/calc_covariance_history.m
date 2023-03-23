function [P_target, P_t] = calc_covariance_history(x, t, stm_t, tcm_time, simparams)




%%%%%%%%%%%%%%%%%%% TODO: ADD DV CALCS TO THIS FUNCTION, AN EASY ADD TO
%%%%%%%%%%%%%%%%%%% EACH LOOP

P_initial = simparams.P_initial;
% maneuverSegments = simparams.maneuverSegments;

% fixed_xfer_duration = simparams.fixed_xfer_duration;

m = simparams.m; % the number of elements per segment
n = simparams.n; % the number of segments

% Reshaping x so each column is a segment initial state and duration
x = reshape(x,m,n);


R = zeros(6,6);
R(4:6,4:6) = simparams.R;
% R = simparams.R; % the TCM execution error covariance

% STM from beginning of trajectory to the state being targeted

if simparams.target_final_maneuver
    final_seg = simparams.maneuverSegments(end);
    target_time = sum(x(7,1:final_seg - 1));
    target_idx = find(target_time == t) ;
    stmN0 = stm_t(:,:,target_idx);
    % Abbreviate t and stm_t to go only until the target (not the end of
    % the trajectory)
    t = t(1:target_idx);
    stm_t = stm_t(:,:,1:target_idx);

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


%% Logic if no TCM (zero length)
if length(tcm_time) == 0
    if nargout == 2
        P_t = tmult(stm_t, tmult(P_i, stm_t, [0 1]));
    end

    P_target = stmN0 * P_i * stmN0';

else % Otherwise, if there are TCMs, do:

    %% Setup for logic choices and STM portions
    % Test which segment has the correction
    
    tcm_idx_logical = logical(sum(t'==tcm_time', 1));
    tcm_t_idx = find(tcm_idx_logical);
    
    % corrSeg = t_s(tcm_idx_logical);
    
    % % Get the nominal state at the correction (DON'T THINK IS NEEDED, BUT
    % SAVING JUST IN CASE)
    % x_tcm = x_t(tcm_idx_logical,:)';
    
    % Extract STMs for calculations
    % A tensor of STMs from the beginning of the trajectory to each correction
    stmC0 = stm_t(:,:,tcm_idx_logical);
    



    
    %% Perform the covariance propagation
    % If the covariance history is being requested (output 2)
    if nargout == 2
        % Compute covariance history
        for i = 1:num_tcm
    
            if i == 1
                P_t(:,:,1:tcm_t_idx(i)) = tmult( stm_t(:,:,1:tcm_t_idx(i)), tmult(P_i, stm_t(:,:,1:tcm_t_idx(i)), [0 1])  );
    
                % Save to use on the next iteration
    %             stmClast0 = stmC0(:,:,i);
    %             stm0Clast = -J * stmClast0' * J;
            else
                %%%%%% RETURN HERE TO VERIFY THE FOLLOWING 3 LINES ARE WORKING
                stmi0_ten = stm_t(:,:,tcm_t_idx(i-1)+1:tcm_t_idx(i));
                stmiClast_ten = tmult(stmi0_ten, stm0C);
                P_t(:,:,tcm_t_idx(i-1)+1:tcm_t_idx(i)) = tmult(stmiClast_ten , tmult(P_t(:,:,tcm_t_idx(i-1)), stmiClast_ten, [0 1])  );
    
    
                % Save to use on the next iteration
    %             stmClast0 = stmC0(:,:,i);
    %             stm0Clast = -J * stmClast0' * J;
            end
    
            
    
            % Perform the TCM update
            stm0C = -J * stmC0(:,:,i)' * J;
            stmNC = stmN0 * stm0C;
    
            T = [-inv( stmNC(1:3,4:6) ) * stmNC(1:3,1:3), -eye(3)];
            N = [zeros(3,6); T];
            IN = eye(6) + N;
    
    
    %         if i == 2
    %             R = zeros(6,6);
    %         end
    
            P_t(:,:,tcm_t_idx(i)) = IN * P_t(:,:,tcm_t_idx(i)) * IN' + R;
    
    
        
        end
    
        % Propagate the remainder of the trajectory dispersion covariance
    
        % need to multiply stmi0 * stm0Clast = stmiClast to propgate the covariance
        % through each time for the remainder of the traj: stmiClast * P_c * stmiClast'
        stmi0_t = stm_t(:,:,tcm_t_idx(i)+1:end);
    %     stm0i_t = tmult(-J, tmult(stmi0_t, J, [1 0]));
    
    %     stmClast0 = stmC0(:,:,end);
    %     stm0Clast = -J * stmClast0' * J;
        stmiClast_t = tmult(stmi0_t, stm0C);
    
        P_t(:,:,tcm_t_idx(i)+1:end) = tmult(  stmiClast_t, tmult(P_t(:,:,tcm_t_idx(i)), stmiClast_t, [0 1])  );
    
    end
    
    R(4:6,4:6) = simparams.R;
    
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
    
        
    end
    
    % Propagate the rest of the way to the target
    % stmClast0 = stmC0(:,:,end);
    % stm0Clast = -J * stmClast0' * J;
    % stmNClast = stmN0 * stm0Clast;
    
    P_target = stmNC * P_i * stmNC';

end


end