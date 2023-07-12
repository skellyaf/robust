function [deltaV_gradient] = calc_deltaV_gradient(x, x_i_f, stm_i, simparams)
%calc_deltaV_gradient Computes and returns the analytical objective function
%gradient
%   Detailed explanation goes here


m = simparams.m;
n = simparams.n;
mu = simparams.mu;
dynSys = simparams.dynSys;

% Reshaping x so each column is a segment initial state and duration
x = reshape(x,m,n);

maneuverSegments = simparams.maneuverSegments;

% Gradient of delta V portion of objective function
deltaV_gradient = zeros(n*m,1);


    
idv = zeros(3,n-1);
    
for i = 1:n
    if i < n
        idv(:,i) = ( x(4:6,i+1) - x_i_f(4:6,i) ) / norm( x(4:6,i+1) - x_i_f(4:6,i) );
    end
    if ismember(i+1, maneuverSegments)
        
        % Partial of deltaV magnitude at end of current segment
        % with respect to initial x vector for segment
        % ex: partial norm ( v0,2 - vf,1 ) / partial x1
        dDVndx_minus = -idv(:,i)' * stm_i(4:6,1:6,i);
        deltaV_gradient(i*m-6:i*m-1) = deltaV_gradient(i*m-6:i*m-1) + dDVndx_minus';

        % Partial of deltaV magnitude at end of current segment
        % with respect to final time for segment 
        % ex: partial norm ( v0,2 - vf,1 ) / partial tf,1
        xdot_xif = stateDot(x_i_f(:,i), mu, dynSys);
        vdot_xif = xdot_xif(4:6);
        deltaV_gradient(i*m,1) = deltaV_gradient(i*m,1) - idv(:,i)' * vdot_xif;
    end

    if ismember(i, maneuverSegments)
        % Partial of deltaV magnitude at beginning of current segment
        % with respect to initial x vector for segment
        % ex: partial norm ( v0,2 - vf,1 ) / partial x2
        deltaV_gradient(i*m-6:i*m-1,1) = deltaV_gradient(i*m-6:i*m-1,1) + ( idv(:,i-1)' * [zeros(3,3), eye(3)] )';
        
        % This one excludes the initial because the deltaV is calculated
        % based on the difference between segment i and i-1

    end
    


end

deltaV_gradient = deltaV_gradient';

end