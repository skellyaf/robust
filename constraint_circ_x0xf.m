function [cin, ceq, cinGrad, ceqGrad] = constraint_circ_x0xf(x, simparams)

% Reshaping x so each column is a segment initial state and duration
m = simparams.m; % the number of elements per segment
n = simparams.n; % the number of segments   
x = reshape(x,m,n);   
nsv = simparams.nsv;
dynSys = simparams.dynSys;

%% Setup
x0 = simparams.x_init;

if ~isfield(simparams,'rdvz_flag')

    x_target = simparams.x_target;
else
    %%%% FLEXIBLE RENDEZVOUS TARGET
    if simparams.rdvz_flag == 1
        total_time = sum(x(7,:));
        simparams.x_target = stateProp(simparams.x0_target, total_time, simparams);
        x_target = simparams.x_target;
    else
        x_target = simparams.x_target;
    end
end
mu = simparams.mu;
maneuverSegments = simparams.maneuverSegments;

fixed_xfer_duration = simparams.fixed_xfer_duration;
if fixed_xfer_duration
    t_f = simparams.tf;
end

if nargout > 2
    outputCGradients = 1;
else
    outputCGradients = 0;
end        


    
    
%% Calculate the dimensions of the equality constraint vector
% Formula is the sum of the following:
%   + 3 + 6: 3 constraints required for a circular orbit of a specific altitude, initial. full state constraint final to a target 
%   + (3 + mod(nsv,6)) * the number of maneuverSegments: only position constraints (and theta_em angle if BR4BP) 
%   + (n+1 -2 (each end constraint) - the # of maneuver segments) * nsv: full state constrained intermediate segments that aren't maneuvers 
%   + 1 if there is a fixed total transfer duration
%   + 2: the initial and final orbits are constrained to be zero length

num_int_full_constraint_nodes = n+1 - 2 - length(maneuverSegments);

% ceq_length = 12 + length(maneuverSegments)*3 + num_int_full_constraint_nodes*6 + logical(simparams.fixed_xfer_duration) + simparams.constrain_flyby_radius;
ceq_length = 3 + 6 + length(maneuverSegments)*(3 + mod(nsv,6)) + num_int_full_constraint_nodes*nsv + logical(simparams.fixed_xfer_duration) + 2;

% Create an empty equality constraint vector
ceq = zeros(1,ceq_length);


%% Calculate the dimensions of the inequality constraint vector

% The only inequality constraints are the dt>=0 (moving forward in time constraint)
% There is one for each segment, so the length is the number of segments

% For the flyby constraint, move the +simparams.constrain_flyby_radius to
% here if want it to be an inequlity constraint, or back to ceq_length if
% want it to be an equality constraint.
% Minus 2 becuase the initial and final segment durations are constrained
% to zero length
% + 1: theta_em_fn constrained to be less than a certain value
cin_length = n - 2 + simparams.constrain_flyby_radius * 2 + 1;
% cin_length = n + simparams.constrain_flyby_radius * 3;
% cin_length = n + simparams.constrain_flyby_radius * 4;  % for debugging

cin = zeros(1,cin_length);

% Preallocate
if outputCGradients
    ceqGrad = zeros(ceq_length, n*m);
    cinGrad = zeros(cin_length, n*m);
end
    
%% Pre-calculate/propagate states/STMs before assembling constraints

if existsAndTrue('target_Pr_constraint_on',simparams)

    traj = createStateStmSttQdQHistory(x(:), simparams);
    [~, deltaVs_nom] = calcDeltaV(x, traj.x_i_f, traj.stm_i, simparams);

    % Assign TCM times to their assigned nodes (ALTERNATIVE METHOD)
    tcm_time = zeros(1, length(simparams.tcm_nodes));
    tcm_idx = zeros(1, length(simparams.tcm_nodes));

    for i = 1:length(simparams.tcm_nodes)                
        tcm_time(i) = sum(x(m,1:simparams.tcm_nodes(i)-1));
    end

    for i = 1:length(tcm_time)
        tcm_idx(i) = find(traj.t == tcm_time(i))';
    end

    [Q_k_km1, dQ_k_km1_dxi, dQ_k_km1_ddti] = calc_Q_events(traj, x, tcm_time, simparams);
    [P_target,~,~,P_k_minus] = calc_covariance_wQ_tcmdv_v3(x, traj, tcm_time, 1, deltaVs_nom, simparams.P_initial, Q_k_km1, simparams);
    [~, ~, ~, dPCkminusdxi, dPCkminusddti] = calc_multiple_tcm_gradient_wQ(x, traj, tcm_time, tcm_idx, P_k_minus, dQ_k_km1_dxi, dQ_k_km1_ddti, deltaVs_nom, simparams);

else
    traj  = createStateStmHistory(x, simparams);

end




%% Now assemble constraints with pre-calculated parameters

neq = 0; % Equality constraint counter
niq = 0; % Inequality constraint counter    

for i = 1:n        
    
    % Extract from segment vector
    x_i_initial = x(1:nsv,i);
    x_i_final = traj.x_i_f(:,i); % used for dynamics constraints        
    delta_t = x(m,i);
    
    % Constraints on the first segment
    if i == 1
%         % Constrain first segment initial position and velocity to
%         % given initial state
% %         ceq(neq+1:neq+nsv) = x_i_initial - x0; % Initial state with theta_em constrained
%         ceq(neq+1:neq+6) = x_i_initial(1:6) - x0(1:6); % Initial state with theta_em UN-constrained
%         
% 
%         if outputCGradients
% %             ceqGrad(neq+1:neq+6,(i-1)*m + 1:i*m) = [eye(nsv), zeros(nsv,1)]; With theta_em constrained 
%             ceqGrad(neq+1:neq+6,(i-1)*m + 1:i*m) = [eye(6), zeros(6,m-6)]; % Without theta_em constrained
%         end
% 
% %         neq = neq + nsv; % Theta_em constrained
%         neq = neq + 6; % With theta_em unconstrained


        % 3 Circular initial orbit constraints of specific radius
        theta_em01 = x_i_initial(nsv);
        xe_0 = 1 - simparams.mub - 1/simparams.a4 * simparams.mu * cos(theta_em01);
        ye_0 = -1/simparams.a4 * simparams.mu *sin(theta_em01);
        r_e = [xe_0; ye_0; 0];

        if outputCGradients
            [ceq(neq+1:neq+3), ceqGrad(neq+1:neq+3,(i-1)*m + 1:i*m)] = circular_initial_orbit_constraints(x_i_initial, r_e, simparams.fixed_initial_radius, simparams.mu_earth_nd, simparams);
        else
            ceq(neq+1:neq+3) = circular_initial_orbit_constraints(x_i_initial, r_e, simparams.fixed_initial_radius, simparams.mu_earth_nd, simparams);
        end
        neq = neq + 3;


    end        
    
    % Constraint to allow an impulsive maneuver at end of first seg and
    % beginning of last seg
    if ismember(i+1,maneuverSegments)

        if i+1 == simparams.n+1 % No final coast segment
            % Only constrain the ith segment final position to the target
            % position
            
            ceq(neq+1:neq+3) = x_i_final(1:3) - x_target(1:3);

            % If BR4BP, no constraint on final theta_em

            if outputCGradients
                % For now, only the gradient wrt the first half of the
                % target position constraint equation (wrt x_i_final)
                ceqGrad(neq+1:neq+3, (i-1)*m + 1:i*m) = [traj.stm_i(1:3,:,i), x_i_final(4:6)];

                % The only relevant partial derivative for x_target is wrt
                % segment durations (it is a fixed propagation from the
                % target's initial position by the total segment duration)

                % Need to loop through and add time sensitivity to the
                % appropriate constraint equation indices, but once we're
                % outside the current for loop. Saving the
                % target_constraint incides below to reference later.

                t_idxs = logical(repmat([zeros(1,nsv), 1], 1, simparams.n));

                ceqGrad(neq+1:neq+3, t_idxs) = ceqGrad(neq+1:neq+3, t_idxs) - ones(3,simparams.n) .* x_target(4:6);


                

            end

            neq = neq + 3;




        else
            % Constrain only position on either side of maneuver nodes for
            % all intermediate nodes    
    
            %  Constrain end position with starting position of next segment
            x_next = x(1:nsv,i+1); % Beginning state of next segment
            ceq(neq+1:neq+3) = x_next(1:3) - x_i_final(1:3);
    
            if outputCGradients
                ceqGrad(neq+1:neq+3, (i-1)*m + 1:i*m) = [-traj.stm_i(1:3,:,i), -x_i_final(4:6)];
                ceqGrad(neq+1:neq+3, i*m + 1:(i+1)*m) = [eye(3), zeros(3,m-3)];
            end    
            neq = neq + 3;

            if strcmp(dynSys, 'br4bp_sb1')
                % Then also constrain the theta_em angle at the end of one
                % segment to the beginning of the next
                ceq(neq+1) = x_next(nsv) - x_i_final(nsv);
                if outputCGradients
                    ceqGrad(neq+1, (i-1)*m + 1:i*m) = [-traj.stm_i(nsv,:,i), -simparams.theta_em_dot];
                    ceqGrad(neq+1, i*m + 1:(i+1)*m) = [zeros(1,nsv-1), 1, 0];
                end   
                neq = neq + 1;
            end
        end
        
    % All other intermediate segments - full state connection constraint
    elseif i < n            
        %  Constrain end position and velocity with beginning of next segment
        x_next = x(1:nsv,i+1); % Beginning state of next segment
        ceq(neq+1:neq+nsv) = x_next(1:nsv) - x_i_final(1:nsv);


        

        if outputCGradients
            ceqGrad(neq+1:neq+nsv, (i-1)*m + 1:i*m) = [-traj.stm_i(:,:,i), -stateDot(x_i_final, simparams)];
            ceqGrad(neq+1:neq+nsv, i*m + 1:(i+1)*m) = [eye(nsv), zeros(nsv,1)];
        end


        neq = neq + nsv;
    elseif i == n % Final segment  
        
%         % Constrain end of final segment to target state
%         ceq(neq+1:neq+nsv) = x_i_final - x_target;
% 
%         
% 
%         if outputCGradients
%             ceqGrad(neq+1:neq+nsv,(i-1)*m + 1:i*m) = [traj.stm_i(:,:,i), stateDot(x_i_final, simparams)];
%         end
% 
%         neq = neq + nsv;

        % 3 Circular target orbit constraint
        theta_emfn = x_i_final(nsv);
        xm_f = 1 - simparams.mub + 1/simparams.a4 * (1-simparams.mu) * cos(theta_emfn);
        ym_f = 1/simparams.a4 * (1-simparams.mu) * sin(theta_emfn);
        
        r_mf = [xm_f; ym_f; 0];

        if outputCGradients
% %             [ceq(neq+1:neq+3), ceqGrad(neq+1:neq+3,(i-1)*m + 1:i*m)] = circular_target_orbit_constraints(x_i_final, r_mf, traj.stm_i(:,:,n), simparams.fixed_final_radius, simparams.mu_moom_nd, simparams);
%             [ceq(neq+1:neq+3), ceqGrad(neq+1:neq+3,(i-1)*m + 1:i*m)] = circular_initial_orbit_constraints(x_i_initial, r_mf, simparams.fixed_final_radius, simparams.mu_moom_nd, simparams);
            
            % circular_lunar_target_constraints
            [ceq(neq+1:neq+6), ceqGrad(neq+1:neq+6,(i-1)*m + 1:i*m)] = circular_lunar_target_constraints(x_i_final, r_mf, traj.stm_i(:,:,n), simparams.fixed_final_radius, simparams.mu_moom_nd, simparams);
        else
% %             ceq(neq+1:neq+3) = circular_target_orbit_constraints(x_i_final, r_mf, traj.stm_i(:,:,n), simparams.fixed_final_radius, simparams.mu_moom_nd, simparams);
%             ceq(neq+1:neq+3) = circular_initial_orbit_constraints(x_i_initial, r_mf, simparams.fixed_final_radius, simparams.mu_moom_nd, simparams);
     
            % circular_lunar_target_constraints
            ceq(neq+1:neq+6) = circular_lunar_target_constraints(x_i_final, r_mf, traj.stm_i(:,:,n), simparams.fixed_final_radius, simparams.mu_moom_nd, simparams);
     
        end


        neq = neq + 6;
                               
    end
    
    % Constraint on time of each segment
    if i == 1 || i == n
        % Initial and final segment duration constrained to zero (only
        % applies to this constraint file with the initial and final
        % circular orbit constraint)
        ceq(neq+1) = delta_t;

        if outputCGradients
            ceqGrad(neq+1, i*m) = 1;
        end
        neq = neq+1;
    else
        % Moving forward in time inequality constraint (each seg time >= 0)
        cin(niq+1) = - delta_t;
    
        if outputCGradients
            cinGrad(niq+1, i*m) = -1;
        end
    
    
        niq = niq + 1;
    end
        
    
                
      
   
    
    
end % end of the main for loop



%% Single constraints

% Constrain sum of delta_t's to equal t_f if fixed_xfer_duration is
% flagged


if fixed_xfer_duration
    total_time = sum(x(m,:));
    ceq(neq+1) = total_time - t_f;
    if outputCGradients
        t_idxs = logical(repmat([zeros(1,nsv), 1], 1, simparams.n));
        ceqGrad(neq+1, t_idxs) = 1;

    end
    neq = neq+1;
end

% Constrain theta_em_0 and theta_em_f between max and min values
%%%%%%%%%%%%%%%%%%%%%%%%%%%% COME BACK HERE
if strcmp(dynSys,'br4bp_sb1')
    % Constrain sin(theta_em_f) to be less than -0.6
    % constraint: sin(theta_em_f) + 0.6 <= 0
    theta_emfn = x_i_final(nsv);
    cin(niq+1) = sin(theta_emfn) + .8;
    if outputCGradients
        cinGrad(niq+1, m*n - 1) = cos(theta_emfn);
    end
    niq = niq + 1;


    % Constrain cos(theta_em_f) to be less than 0.25
    % constraint: cos(theta_em_f) - 0.25 <= 0
%     cin(niq+1) = cos(theta_emfn) - 0.1;
    cin(niq+1) = cos(theta_emfn) - 0;
    if outputCGradients
        cinGrad(niq+1, m*n - 1) = -sin(theta_emfn);
    end
    niq = niq+1;
end


% Option 3: Find perilune, constrain the distance from the moon to be
% greater than or equal to the limit
if simparams.constrain_flyby_radius

    % Find the closest approach to the moon (perilune)
    %%%%% MODIFIED TO WORK WITH ONLY BR4BP
    r_sc_t = traj.x_t(:,1:3)';
    theta_t = traj.x_t(:,7);
    % Find the position of the earth as a function of time and theta
    xe_t = 1 - simparams.mub - 1/simparams.a4 * simparams.mu * cos(theta_t);
    ye_t = -1/simparams.a4 * simparams.mu *sin(theta_t);
    re_t = [xe_t, ye_t, zeros(length(xe_t),1)]';
    r_sc_e_t = r_sc_t - re_t;
    d_e_sc_t = vecnorm(r_sc_e_t);
    [d_perilune, idx_perilune] = min(d_e_sc_t);

    % Constrain the closest approach to be > simparams.flyby_radius (simparams.flyby_radius - d_perilune < 0) 
    cin(niq+1) = simparams.flyby_radius - d_perilune;

    % Perilune distance inequality constraint gradient
    if outputCGradients
        % Find the segment
        seg_perilune = traj.t_s(idx_perilune);

        %
        i_m_sc = r_sc_e_t(:,idx_perilune) / d_perilune;

        % Find the index at the beginning of the segment with perilune
%         pSeg_indices = t_s == seg_perilune
%         idx_start_pSeg = 



        t0_pSeg = sum(x(m,1:seg_perilune-1));
        idx_pSeg = find(traj.t==t0_pSeg);


        stm_perilune_seg0 = dynCellCombine(traj.t, traj.t_s, idx_pSeg, idx_perilune, simparams, traj.stm_t_i);
        k = seg_perilune;
        cinGrad(niq+1, (k-1)*m + 1 : (k-1)*m + nsv) = -i_m_sc' * stm_perilune_seg0(1:3,:);







    end



    % Increment the number of inequality constraints
    niq = niq+1;

end

% Target position covariance RSS inequality constraint
if existsAndTrue('target_Pr_constraint_on',simparams)

    [event_times, event_indicator] = define_events_v2(x(:), traj.t, tcm_time, simparams);

%     event_idx_logical = logical(sum(traj.t'==event_times', 1));    
%     event_idxs = find(event_idx_logical);

    M = [eye(3), zeros(3,nsv-3)];

    k_last = simparams.start_P_growth_node;


    for q = 1:length(simparams.P_constrained_nodes)
        % Get the event index for the covariance at (minus) the P
        % constrained node
        t_k = sum(x(m,1:simparams.P_constrained_nodes(q)-1));
        [~, k] = find(event_times == t_k);


        % Ensuring it is a plain nominal maneuver (0) or a corrected
        % nominal maneuver (3) w/the following assert:
        assert(event_indicator(k) == 0 || event_indicator(k) == 3);

        % Get the covariance of interest (minus):
        P_k = P_k_minus(:,:,k);

        % Create position dispersion constraint equation there 
        % CURRENTLY THE SAME MAX P_R AT EACH K...COULD BE AN
        % UPGRADE/CUSTOMIZATION
        cin(niq+1) = sqrt(trace(M*P_k*M')) - simparams.P_max_r;

        if outputCGradients    
            % Analytical Pr constraint gradient
            for i = 1:n
                for j = 1:nsv
                    cinGrad(niq+1, (i-1)*m + j) = 1/2 * trace(M*P_k*M')^(-1/2) * trace(M*dPCkminusdxi(:,:,j,k,i)*M');
                end
                    cinGrad(niq+1, i*m) = 1/2 * trace(M*P_k*M')^(-1/2) * trace(M*dPCkminusddti(:,:,k,i)*M');
            end
        end



        niq = niq + 1; % Increment inequality constraints total


    end






    
% % % %     cin(niq + 1) = sqrt(trace(M*P_target*M')) - simparams.P_max_r;
% % % % 
% % % %     if outputCGradients
% % % %         [event_times] = define_events_v2(x(:), traj.t, tcm_time, simparams);
% % % %         m = length(event_times); % Total number of TCM modifying events / POTENTIALLY UNNECESSARY...
% % % %         % Adding an assert here that the third dimension of P_k_minus
% % % %         % equals m. If it never throws an error, can just index 'end' on
% % % %         % the third for P_k_minus instead of m.
% % % % 
% % % % 
% % % %         % Analytical Pr constraint gradient
% % % %         for i = 1:simparams.n
% % % %             for j = 1:6
% % % % %                 dPr_an((i-1)*7 + j) = 1/2 * trace(M*P_target*M')^(-1/2) * -trace(M*dPCkminusdxi(:,:,j,m,i)*M');
% % % %                 cinGrad(niq+1, (i-1)*7 + j) = 1/2 * trace(M*P_target*M')^(-1/2) * trace(M*dPCkminusdxi(:,:,j,m,i)*M');
% % % %             end
% % % % %                 dPr_an(i*7) = 1/2 * trace(M*P_target*M')^(-1/2) * -trace(M*dPCkminusddti(:,:,m,i)*M');
% % % %                 cinGrad(niq+1, i*7) = 1/2 * trace(M*P_target*M')^(-1/2) * trace(M*dPCkminusddti(:,:,m,i)*M');
% % % %         end
% % % % 
% % % % 
% % % % 
% % % % 
% % % %     end
% COMMENTED THE LINE BELOW - IT WAS UNCOMMENTED BEFORE, WHICH I BELIEVE WAS
% A BUG:
% % % %     niq = niq + 1; % Increment inequality constraints total
end



%% 
if outputCGradients
    cinGrad = cinGrad';
    ceqGrad = ceqGrad';
end

end

