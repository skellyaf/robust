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

% Create an empty equality constraint vector
ceq = zeros(1,ceq_length);

if outputCGradients
    ceqGrad = zeros(ceq_length, n*m);
end

%% Calculate the dimensions of the inequality constraint vector

% The only inequality constraints are the dt>=0 (moving forward in time constraint)
% There is one for each segment, so the length is the number of segments
cin_length = n;

cin = zeros(1,cin_length);

if outputCGradients
    cinGrad = zeros(cin_length, n*m);
end
    
%% Pre-calculate/propagate states/STMs before assembling constraints

% Create the state and stm history:
[stm_i, x_i_f, ~, stm_t, t] = createStateStmHistory(x, simparams);

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
    r_b = [1+simparams.mu; 0; 0]; % Position of the moon
    r_n = x(1:3,simparams.flyby_node); % Position of the constrained node
    r_d = r_n - r_b;
    d = vecnorm(r_d);
    ceq(neq+1) = d - simparams.flyby_radius;
    

    % Gradient addition
    if outputCGradients
        i_d = r_d / d;
        k = simparams.flyby_node;
        ceqGrad(neq+1, (k-1)*7 + 1 : (k-1)*7 + 3) = i_d'; % Vector magnitude partial derivative
    end

    neq = neq+1;
end







%% 
if outputCGradients
    cinGrad = cinGrad';
    ceqGrad = ceqGrad';
end

end

