function [deltaV, deltaVs, deltaV_gradient] = calcDeltaV(x, x_i_f, stm_i, simparams)
%calcDeltaV Calculates impulsive delta V between segment final velocities
%and the initial velocities of the next segment
% maneuverSegments identify the segment number with a maneuver at the
% beginning of it




m = simparams.m; % the number of elements per segment
n = simparams.n; % the number of segments
maneuverSegments = simparams.maneuverSegments;
mu = simparams.mu;
dynSys = simparams.dynSys;

% Reshaping x so each column is a segment initial state and duration
x = reshape(x,m,n);


vi = zeros(3,length(maneuverSegments));
vf = zeros(3,length(maneuverSegments));

if nargout > 2
    deltaV_gradient = zeros(m, n);
end

% performing outside this function
% if isfield(simparams,'rdvz_flag')
%     % The target state is flexible and depends on each segment duration
%     %%%% FLEXIBLE RENDEZVOUS TARGET
%     total_time = sum(x(7,:));
%     simparams.x_target = stateProp(simparams.x0_target, total_time, simparams);
% end





for j = 1:length(maneuverSegments)



    if maneuverSegments(j) == 1

        vi(:,1) = x(4:6, 1);
        vf(:,1) = simparams.x_init(4:6);

        % Gradient calc - delta V at the first node / beginning of 1st seg
        if nargout > 2
            i = maneuverSegments(j);
            dv = vi(:,j) - vf(:,j);
            idv = dv / norm(dv);
            dDVndx = idv' * [zeros(3,3), eye(3)];
            deltaV_gradient(1:6, i) = dDVndx;
        end

    elseif maneuverSegments(j) == n+1

        vi(:,j) = simparams.x_target(4:6,:);
        vf(:,j) = x_i_f(4:6,n);

        % Gradient calc - delta V occurs at end of final segment
        if nargout > 2
            i = maneuverSegments(j);
            dv = vi(:,j) - vf(:,j);
            idv = dv / norm(dv);
            
            dDVndx = - idv' * stm_i(4:6,:,n);
            deltaV_gradient(1:6,n) = dDVndx;

            xdot_xif = stateDot(x_i_f(:,n), mu, dynSys);
            deltaV_gradient(7,n) = -idv' * xdot_xif(4:6);

        end

    else

        vi(:,j) = x(4:6, maneuverSegments(j));
        vf(:,j) = x_i_f(4:6, maneuverSegments(j)-1);

        % Gradient calc - delta V at any intermediate segments
        if nargout > 2
            i = maneuverSegments(j);
            dv = vi(:,j) - vf(:,j);
            idv = dv / norm(dv);

            dDVndx_minus = -idv' * stm_i(4:6,:,i-1);
            deltaV_gradient(1:6, i-1) = dDVndx_minus;

            xdot_xif = stateDot(x_i_f(:,i-1), mu, dynSys);
            dDVnddt_minus = -idv' * xdot_xif(4:6);
            deltaV_gradient(7, i-1) = dDVnddt_minus;
            
            dDVndx = idv' * [zeros(3,3), eye(3)];
            deltaV_gradient(1:6, i) = dDVndx;
        end

    end
end



if nargout > 2
    if isfield(simparams,'rdvz_flag')
        if simparams.rdvz_flag == 1
            % The target state is flexible and depends on each segment duration
            %%%% FLEXIBLE RENDEZVOUS TARGET
    %         total_time = sum(x(7,:));
    %         simparams.x_target = stateProp(simparams.x0_target, total_time, simparams);
            x_dot_target = stateDot(simparams.x_target, mu, dynSys);
    
            dv = vi(:,end) - vf(:,end);
            idv = dv / norm(dv);
    
            deltaV_gradient(7,:) = deltaV_gradient(7,:) + idv' * x_dot_target(4:6) * ones(1,n);

        end


    end
    deltaV_gradient = deltaV_gradient(:)';

end









% if ismember(0, maneuverSegments-1) && ismember(n+1, maneuverSegments)
%     assert(0,'need more functionality')
% 
% elseif ismember(0, maneuverSegments-1)
%     vi = x(4:6,maneuverSegments);
%     msm1 = maneuverSegments-1;
%     vf = [simparams.x_init(4:6), x_i_f(4:6,msm1(2:end)) ];
% 
% elseif ismember(n+1, maneuverSegments)
%     if isfield(simparams,'rdvz_flag')
%         vi = [x(4:6,maneuverSegments(1:end-1)), x_i_f(4:6,end)];
%         vf = [x_i_f(4:6, maneuverSegments(1:end-1)-1), simparams.x_target(4:6,:)];
%         
%     else
%         assert(0,'ERROR, SHOULDNT BE HERE OR NEED FUNCTIONALITY')
%     end
% else
%     vi = x(4:6,maneuverSegments);
%     vf = x_i_f(4:6,maneuverSegments-1);
%     
% end

deltaVs = vi-vf;
deltaV = sum(vecnorm(deltaVs));

end