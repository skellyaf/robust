function [stm_i, x_i_f, x_t, stm_t, t] = createStateStmHistory(x, simparams)
%createStateStmHistory Returns the history of the state and STM via numeric
%integration of dynamics

mu = simparams.mu;

% The correction is not part of the multi-segment trajectory vector x (as
% it has been in some instances)

m = simparams.m; % the number of elements per segment
n = simparams.n; % the number of segments

% Reshaping x so each column is a segment initial state and duration
x = reshape(x,m,n);
%% Propagate entire trajectory, save dynamics at each time step

for i = 1:n

    x_i_initial = x(1:6,i);
    delta_t = x(7,i);

    %%%%%%%%%%%% WHAT IF DELTA_T IS NEGATIVE? HOW TO ALLOW THAT?  REQUIRES
    %%%%%%%%%%%% THOUGHT (one example resulted in convergence to an
    %%%%%%%%%%%% infeasible point)
    if delta_t ~= 0
        [x_i_final, stm_i(:,:,i), xstm_t_i, t_i] = stateStmProp(x_i_initial, delta_t, simparams);
        x_i_f(:,i) = x_i_final;

        if i == 1
            % State history
            x_t = xstm_t_i(:,1:6);
            % STM history
            stm_t = reshape( xstm_t_i(:,7:42)',6,6,[] );
            % Time history
            t = t_i;


        else
            % Append to history structure
            % Time (exclude first time element, it duplicates the final of
            % the previous)
            t = [t; t_i(2:end) + sum(x(7,1:i-1))];
            % State
            x_t = [x_t; xstm_t_i(2:end,1:6)];

            % Reshape new STM history tensor
            stm_t_i = reshape( xstm_t_i(:,7:42)',6,6,[] );

            % Combine with previous final STM to continue history from
            % beginning of trajectory
            
            if isempty(stm_t)
                stm_t = stm_t_i;
            else

                stm_t(:,:,end+1:end+size(stm_t_i,3)-1) = tmult(stm_t_i(:,:,2:end),stm_t(:,:,end));      
            end


        end

    else
        x_i_final = x_i_initial;
        x_i_f(:,i) = x_i_final;
        stm_i(:,:,i) = eye(6);

        if i == 1
            
            stm_t = eye(6);
            x_t = [x_i_initial'];
            t = 0;
        end
    end

end







end