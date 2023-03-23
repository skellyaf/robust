function [t_c_opt, min3sigmaDVC, corrGrad, dvC3sigma, t, t_zerox, stm_t, testVals, x_t, min_3sigma_dvf] = optimal_single_tcm(x, P_initial, simparams, varargin)
%optimal_single_tcm calculates the optimal time to perform a single TCM
%along a nominal trajectory
%   By calculating the gradient of the magnitude of the TCM with respect to
%   the time it is performed, an extrema exists at zero crossings. This
%   function calculates this partial derivative by performing a single
%   propagation from the beginning to end of the trajectory and saving a
%   time history of the states and STMs (at the same interval that the
%   variable step integrator selects for numerical integration). In this
%   way, the trajectory dynamics are applied quickly versus performing a
%   numerical integration for each test case.
%
%   Inputs are the trajectory parameter vector (x), the initial state
%   dispersion covariance matrix (P_initial), and the system gravitational
%   parameter (mu)
%
%   Outputs are the optimal time to perform the correction (t_c_opt) and
%   probably more as required/developed. ADD HERE IF REQUIRED

%% Parse segments
% # of parameters per segment
mu = simparams.mu;
m = 7;
if ~isempty(varargin)
    stmOutIdx = varargin{1};
end

% Calculate the number of corrections
% This works currently as long as there aren't m or more corrections
num_corrections = mod(length(x),m);
n = ( length(x) - num_corrections ) / m;

% Array of corrective maneuver times that is num_corrections long
t_corr = x(end-num_corrections+1:end);

% Reshaping x so each column is a segment initial state and duration
x = reshape(x(1:end-num_corrections),m,n);

testVals = []; % a structure for storing test values along trajectory to try and find out what other indicators correspond to the min TCM


%% Propagate entire trajectory, save dynamics at each time step

for i = 1:n

    x_i_initial = x(1:6,i);
    delta_t = x(7,i);

    if delta_t > 0
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
            % Time
            t = [t; t_i + sum(x(7,1:i-1))];
            % State
            x_t = [x_t; xstm_t_i(:,1:6)];

            % Reshape new STM history tensor
            stm_t_i = reshape( xstm_t_i(:,7:42)',6,6,[] );

            % Combine with previous final STM to continue history from
            % beginning of trajectory

            stm_t(:,:,end+1:end+size(stm_t_i,3)) = tmult(stm_t_i,stm_t(:,:,end));          

        end

    else
        x_i_final = x_i_initial;
        x_i_f(:,i) = x_i_final;
        stm_i(:,:,i) = eye(6);
    end

end

% Making sure history sizes match
assert(size(t,1)==size(stm_t,3));

%% Calculate time gradient
corrGrad = zeros(size(t,1),1);
dvC3sigma = zeros(size(t,1),1);
dv_vf_3sigma = zeros(size(t,1),1);

% Start sweep at zero
t_c = 0;

% Initialize stm_C0 to be identity (starts at zero)
stmC0 = eye(6);

% Propagate entire trajectory once to get stmNC
stmNC = eye(6);

stmN0 = stm_t(:,:,end);

% Symplectic unit matrix
J = [zeros(3,3), eye(3,3); -eye(3,3), zeros(3,3)];

for i = 1:size(t,1)

    t_c = t(i);
    x_c = x_t(i,:)';
    stmC0 = stm_t(:,:,i);
    stm0C = -J * stmC0' * J; % symplectic inverse
    stmNC = stmN0 * stm0C;

    T = [-inv( stmNC(1:3,4:6) ) * stmNC(1:3,1:3), -eye(3)];


    % Partial of correction wrt time only calculation
    dstmNCdtc = -stmNC * r2bp_A_matrix(x_c, mu);
    dTdtc = T_partial(stmNC,dstmNCdtc);
    dstmC0dtc = r2bp_A_matrix(x_c, mu) * stmC0;
    corrGrad(i) = sigDvPartial(T,stmC0,P_initial,dTdtc,dstmC0dtc);   

    % Magnitude of correction
    dvC3sigma(i) = 3 * sqrt( trace( T * stmC0 * P_initial * stmC0' * T' ) );

    % Magnitude of correction of final velocity dispersion
    W = [stmNC(4:6,4:6), -stmNC(4:6,1:3)];
    L = [W*T', zeros(3,3)];
    dv_vf_3sigma(i) = 3 * sqrt( trace( L * stmC0 * P_initial * stmC0' * L' ) );

    % Testing other test vals
    %%%% STILL THINKING ABOUT WHAT TO TEST

% % % % %     eigstmNC = eig(stmNC);
% % % % % 
% % % % %     for evidx = 1:6
% % % % %         testVals(evidx,i) = eigstmNC(ev);
% % % % %     end
% % % % % 
% % % % %     if find(stmOutIdx==i)
% % % % %         testVals(end+1,1:36) = stmC0(:);
% % % % %     end



    

end






zci = @(v) find(diff(sign(v)));

ZeroX = @(x0,y0,x1,y1) x0 - (y0.*(x0 - x1))./(y0 - y1);

% Compute the zero crossings. 
% Also include the first index as it may be optimal for the interval.

zxidx = [1; zci(corrGrad)];

t_zerox = [];
corrMag1 = [];
% corrMag2 = [];


for m = 1:length(zxidx)
    idx1 = zxidx(m);
    idx2 = idx1+1;
    t_1 = t(idx1);
    t_2 = t(idx2);

    y1 = corrGrad(idx1);
    y2 = corrGrad(idx2);

    % Are the slopes the same sign or different signs before, across, and after the
    % correction? (ie is the sign change a zero crossing or asymptote?)

    slope2 = corrGrad(idx2) - corrGrad(idx1);
    if idx1 ~= 1
        slope1 = corrGrad(idx1) - corrGrad(idx1-1);
    else
        slope1 = 1;
    end
    
    

    if slope1 == 0
        slope2 = slope1;
    end

    if idx2+1 < length(corrGrad)
        slope3 = corrGrad(idx2+1) - corrGrad(idx2);
        if slope3 == 0
            slope3 = slope2;
        end
        
    else 
        slope3 = slope2; % Go ahead and evaluate the point
    end

    % If it is zero crossing or the first point
    if abs(  sign(slope1) + sign(slope2) + sign(slope3)   ) == 3 || idx1 == 1
        if idx1 == 1
            t_zerox(end+1) = t(idx1);
            stmC1 = eye(6);
        else
            t_zerox(end+1) = ZeroX(t_1,y1,t_2,y2);
            % Propagate STM from t_1 to t_zerox
            [~, stmC1] = stateStmProp(x_i_initial, t_zerox(end) - t_1, simparams, eye(6));
        end

        
        

        stmC0 = stmC1 * stm_t(:,:,idx1);
        stm0C = -J * stmC0' * J; % symplectic inverse
        stmNC = stmN0 * stm0C;

        T = [-inv( stmNC(1:3,4:6) ) * stmNC(1:3,1:3), -eye(3)];


        % Calculate correction DV at t_zerox
        corrMag1(end+1) = 3 * sqrt( trace( T * stmC0 * P_initial * stmC0' * T' ) );


        % Calculate correction magnitude
%         corrMag = 




%         x_currentCorrection = [x(:); t_zerox(end)];
%         [corrMag2(end+1)] = xfer_obj_corr_grad(x_currentCorrection,P_initial, mu);



    end


    
end

[min3sigmaDVC, minIdx] = min(corrMag1);

min_3sigma_dvf = min(dv_vf_3sigma);

t_c_opt = t_zerox(minIdx);




end