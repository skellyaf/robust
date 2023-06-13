function [cin, ceq, cinGrad, ceqGrad] = constraint_min_tcm(x, simparams)

%% Setup
    x0 = simparams.x_init;
    
    x_target = simparams.x_target;
    mu = simparams.mu;
    maneuverSegments = simparams.maneuverSegments;

    fixed_xfer_duration = simparams.fixed_xfer_duration;
    if fixed_xfer_duration
        t_f = simparams.tf;
    end

    if nargout == 4
        outputCGradients = 1;
    else
        outputCGradients = 0;
    end        
    
    m = simparams.m; % the number of elements per segment
    n = simparams.n; % the number of segments   
    
    
%% Calculate the dimensions of the equality constraint vector
% Formula is the sum of the following:
%   + 12: always constrained position and velocity for initial and final state 
%   + 3 * the number of maneuverSegments: only position constraints
%   + (n+1 -2 (each end constraint) - the # of maneuver segments) * 6: constrained intermediate segments that aren't maneuvers 
%   + 1 if there is a fixed total transfer duration

num_int_full_constraint_nodes = n+1 - 2 - length(maneuverSegments);

ceq_length = 12 + length(maneuverSegments)*3 + num_int_full_constraint_nodes*6 + logical(simparams.fixed_xfer_duration) + simparams.constrain_flyby_radius;
% ceq_length = 12 + length(maneuverSegments)*3 + num_int_full_constraint_nodes*6 + logical(simparams.fixed_xfer_duration);

% Create an empty equality constraint vector
ceq = zeros(1,ceq_length);

if outputCGradients
    ceqGrad = zeros(ceq_length, n*m);
end

%% Calculate the dimensions of the inequality constraint vector

% The only inequality constraints are the dt>=0 (moving forward in time constraint)
% There is one for each segment, so the length is the number of segments

% For the flyby constraint, move the +simparams.constrain_flyby_radius to
% here if want it to be an inequlity constraint, or back to ceq_length if
% want it to be an equality constraint.
cin_length = n + simparams.constrain_flyby_radius;
% cin_length = n + simparams.constrain_flyby_radius * 3;
% cin_length = n + simparams.constrain_flyby_radius * 4;  % for debugging

cin = zeros(1,cin_length);

if outputCGradients
    cinGrad = zeros(cin_length, n*m);
end
    
%% Pre-calculate/propagate states/STMs before assembling constraints

% Create the state and stm history:
[stm_i, x_i_f, x_t, ~, t, t_s, stm_t_i] = createStateStmHistory(x, simparams);

% Reshaping x so each column is a segment initial state and duration
x = reshape(x,m,n);   


%% Now assemble constraints with pre-calculated parameters

neq = 0; % Equality constraint counter
niq = 0; % Inequality constraint counter    

for i = 1:n        
    
    % Extract from segment vector
    x_i_initial = x(1:6,i);
    x_i_final = x_i_f(:,i); % used for dynamics constraints        
    delta_t = x(7,i);
    
    % Constraints on the first segment
    if i == 1
        % Constrain first segment initial position and velocity to
        % given initial state
        ceq(neq+1:neq+6) = x_i_initial - x0;
        

        if outputCGradients
            ceqGrad(neq+1:neq+6,(i-1)*7 + 1:i*7) = [eye(6), zeros(6,1)];
        end

        neq = neq + 6;
    end        
    
    % Constraint to allow an impulsive maneuver at end of first seg and
    % beginning of last seg
    if ismember(i+1,maneuverSegments)
        %  Constrain end position with starting position of next segment
        x_next = x(1:6,i+1); % Beginning state of next segment
        ceq(neq+1:neq+3) = x_next(1:3) - x_i_final(1:3);

        if outputCGradients
            ceqGrad(neq+1:neq+3, (i-1)*7 + 1:i*7) = [-stm_i(1:3,:,i), -x_i_final(4:6)];
            ceqGrad(neq+1:neq+3, i*7 + 1:(i+1)*7) = [eye(3), zeros(3,4)];
        end



        neq = neq + 3;
    
    %%%%%%%%%%%%%% ONE TIME TESTING CONSTRAINT - FORCE ALL PLANE CHANGE
    %%%%%%%%%%%%%% TO OCCUR ON FIRST BURN FOR SCENARIO 2
        if simparams.constrain_dv1_inclination_change    
            %%%%% try 3 - enforce r_tgt x v_tgt to be coplanar with r_xfer
            %%%%% x v_xfer. i.e., create unit vectors for each plane, dot
            %%%%% the unit vectors should be equal to 1
            if ismember(i+1,maneuverSegments(end)) % this if statement just makes sure it only creates one constraint per round
                rf = simparams.x_target(1:3);
                vf = simparams.x_target(4:6);
                rfxvf = cross(rf,vf);
                i_tgt = rfxvf/norm(rfxvf);

                r2f = x_i_f(1:3,2);
                v2f = x_i_f(4:6,2);
                r2fxv2f = cross(r2f,v2f);
                i_xfer = r2fxv2f / norm(r2fxv2f);

                ceq(neq+1) = dot( i_tgt, i_xfer ) - 1;

                % NO ANALYTICAL GRADIENT WRITTEN FOR THIS CONSTRAINT
                neq = neq+1;        
            end    
        end
        
    % All other intermediate segments - full state connection constraint
    elseif i < n            
        %  Constrain end position and velocity with beginning of next segment
        x_next = x(1:6,i+1); % Beginning state of next segment
        ceq(neq+1:neq+6) = x_next(1:6) - x_i_final(1:6);


        

        if outputCGradients
            ceqGrad(neq+1:neq+6, (i-1)*7 + 1:i*7) = [-stm_i(:,:,i), -stateDot(x_i_final, mu, simparams.dynSys)];
            ceqGrad(neq+1:neq+6, i*7 + 1:(i+1)*7) = [eye(6), zeros(6,1)];
        end


        neq = neq + 6;
    end
    
    % Moving forward in time inequality constraint (each seg time >= 0)
    cin(niq+1) = - delta_t;

    if outputCGradients
        cinGrad(i, i*m) = -1;
    end


    niq = niq + 1;
    

            
    % Final segment - constrained end to match target
    if i == n
        % Constrain end of final segment to target state
        ceq(neq+1:neq+6) = x_i_final - x_target;

        if outputCGradients
            ceqGrad(neq+1:neq+6,(i-1)*7 + 1:i*7) = [stm_i(:,:,i), stateDot(x_i_final, mu, simparams.dynSys)];
        end

        neq = neq+6;
                               
    end  
   
    
    
end % end of the main for loop

%% Single constraints

% Constrain sum of delta_t's to equal t_f if fixed_xfer_duration is
% flagged

total_time = sum(x(7,:));
if fixed_xfer_duration
    ceq(neq+1) = total_time - t_f;
    neq = neq+1;
end

% Powered flyby node distance from moon constraint

if simparams.constrain_flyby_radius
    r_b = [1-simparams.mu; 0; 0]; % Position of the moon
    r_n = x(1:3,simparams.flyby_node); % Position of the constrained node
    r_d = r_n - r_b;
    d = vecnorm(r_d);
    % As an equality constraint
%     ceq(neq+1) = d - simparams.flyby_radius;
    % As an inequality constraint
    cin(niq+1) = simparams.flyby_radius - d;
    

    % Gradient addition
    if outputCGradients
        i_d = r_d / d;
        k = simparams.flyby_node;
        % Gradient if it is an equality constraint
%         ceqGrad(neq+1, (k-1)*7 + 1 : (k-1)*7 + 3) = i_d'; % Vector magnitude partial derivative

        % Gradient if it is an inequality constraint
        cinGrad(niq+1, (k-1)*7 + 1 : (k-1)*7 + 3) = - i_d';
    end

%     neq = neq+1;
    niq = niq+1;
end

% Powered flyby node constraint options to ensure other points are not passing within the lunar surface
% Option 1: the vector from the moon to the spacecraft (r_m_sc) is
% orthogonal to the velocity vector (with or without the delta V?)
%%%%%% not really working with the delta V added in...trying subtracting it

if simparams.constrain_flyby_radius

    x_flyby = x(1:6,simparams.flyby_node); % State at powered flyby
    r_b = [1-simparams.mu; 0; 0]; % Position of the moon

    xf_pre_flyby = x_i_f(1:6,simparams.flyby_node - 1); % State at the node prior to the flyby node
    stm_prev = stm_i(:,:,simparams.flyby_node - 1);








%     dv_to_subtract = x_flyby(4:6) - x_i_f(4:6,simparams.flyby_node-1)

    [apse_constraint_eqn, apse_constraint_gradient] = apse_constraint_prev_vel_S(x_flyby, xf_pre_flyby, r_b, stm_prev, simparams);
%     [apse_constraint_eqn, apse_constraint_gradient] = apse_constraint(x_flyby, r_b);

    % Apse constraint as an equality constraint
    ceq(neq+1) = apse_constraint_eqn;

    if outputCGradients
        k = simparams.flyby_node;
    
        % Apse equality constraint gradient
        ceqGrad(neq+1,(k-2)*7 + 1 : k*7) = apse_constraint_gradient;
    end

    neq = neq + 1;

end





% Option 2: enforce that the position at the powered flyby is closer to the
% moon than the state (not node) immediately before it or after it
%%%% having trouble getting the option 2 numerical gradients to match...

% if simparams.constrain_flyby_radius
% 
%     % Find the time of the powered flyby
%     t_flyby = sum(x(7,1:simparams.flyby_node-1));
%     idx_flyby = find(t==t_flyby);
% 
%     rm = [1-mu; 0; 0];
% 
%     r_sc_flyby = x(1:3,simparams.flyby_node);
% 
%     r_m_sc_flyby =  r_sc_flyby - rm;
%     
%     dist_flyby = vecnorm(r_m_sc_flyby);
% 
%     % How many indices are in the segment after the flyby
%     num_ind_before = sum(simparams.flyby_node-1 == t_s);
%     num_ind_before = 25;
% %     idx_before = floor(idx_flyby - .2 * num_ind_before); % Go 1/4 of the way into the segment
%     idx_before = floor(idx_flyby - num_ind_before); % testing a gradient bug
% 
%     r_sc_beforeFlyby = x_t(idx_before, 1:3)';
%     r_m_sc_beforeFlyby = r_sc_beforeFlyby - rm;
%     dist_beforeFlyby = vecnorm(r_m_sc_beforeFlyby);
% 
% 
%     % How many indices are in the segment after the flyby
%     num_ind_after = sum(simparams.flyby_node == t_s);
%     idx_after = floor(idx_flyby + .2 * num_ind_after); % Go 1/4 of the way into the segment
% 
% 
%     r_sc_afterFlyby = x_t(idx_after, 1:3)';
% 
%     r_m_sc_afterFlyby = r_sc_afterFlyby - rm;
%     dist_afterFlyby = vecnorm(r_m_sc_afterFlyby);
% 
% 
% 
% 
%     % Distance < distance_before
%     cin(niq+1) = dist_flyby - dist_beforeFlyby;
% 
%     % Distance < distance_after constraint
%     cin(niq+2) = dist_flyby - dist_afterFlyby;
% 
% 
% 
%     %%%% Debugging - 3 constraints to test which is messing up
% 
% 
% 
% %     cin(niq+1) = dist_beforeFlyby;
% %     cin(niq+2) = dist_afterFlyby;
% %     cin(niq+3) = dist_flyby;
% 
% 
%     %%%% Debugging
%     
% 
%     
% 
% 
% 
% 
%     % Option 2 gradient
%     if outputCGradients
%         k = simparams.flyby_node;
%         % Find the index for the first node of the segment before the flyby segment
%         idx_pre_flyby = find(t_s==k-1);
%         idx_pre_flyby_start = idx_pre_flyby(1);
% 
%         %%%%%% WEIRD THING HAPPENED HERE...NEED TO REMEMBER IF THAT IS HOW
%         %%%%%% I SET IT UP ON PURPOSE. WHEN EVALAUTING t_s(t==sum(x(7,1:simparams.flyby_node-2)))
%         %%%%%% , THE TIME ALONG THE TRAJECTORY AT THE END(?) OF SEGMENT 12,
%         %%%%%% IT IS SLIGHTLY DIFFERENT THAN THE TIME AT THE BEGINNING OF
%         %%%%%% SEGMENT 13: t(idx_pre_flyby_start)
% 
%         
%         i_r_m_sc_flyby = r_m_sc_flyby / dist_flyby;
%         i_r_m_sc_afterFlyby = r_m_sc_afterFlyby / dist_afterFlyby;
%         i_r_m_sc_beforeFlyby = r_m_sc_beforeFlyby / dist_beforeFlyby;
%         
% %         stmCkClast = dynCellCombine(t, t_s, tClast_idx, tCk_idx, simparams, stm_t_i);
% 
%         stm_after_flyby = dynCellCombine(t, t_s, idx_flyby, idx_after, simparams, stm_t_i);
%         stm_flyby_before = dynCellCombine(t, t_s, idx_pre_flyby_start, idx_before, simparams, stm_t_i);
% 
%         % Partial derivative wrt previous segment
%         cinGrad(niq+1, (k-2)*7 + 1 : (k-2)*7 + 6) =  i_r_m_sc_beforeFlyby' * stm_flyby_before(1:3,:);
% 
%         % Partial derivative wrt current segment
%         cinGrad(niq+1, (k-1)*7 + 1 : (k-1)*7 + 6) = i_r_m_sc_flyby' * [eye(3), zeros(3,3)];
%         cinGrad(niq+2, (k-1)*7 + 1 : (k-1)*7 + 6) = i_r_m_sc_flyby' * [eye(3), zeros(3,3)] + i_r_m_sc_afterFlyby' * stm_after_flyby(1:3,:);
% 
% 
% 
% 
% 
% %         % d dist_beforeflyby / d prev segment
% %         cinGrad(niq+1, (k-2)*7 + 1 : (k-2)*7 + 6) = i_r_m_sc_beforeFlyby' * stm_flyby_before(1:3,:);
% % 
% %         % d dist_afterFlyby / d curr segment
% %         cinGrad(niq+2, (k-1)*7 + 1 : (k-1)*7 + 6) = i_r_m_sc_afterFlyby' * stm_after_flyby(1:3,:);
% % 
% %         % d dist_flyby / d curr segment
% %         cinGrad(niq+3, (k-1)*7 + 1 : (k-1)*7 + 6) = i_r_m_sc_flyby' * [eye(3), zeros(3,3)];
%         
%     end
% 
% 
% 
% 
%     niq = niq + 2;
% end








%% 
if outputCGradients
    cinGrad = cinGrad';
    ceqGrad = ceqGrad';
end

end

