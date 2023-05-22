function [ T, X_return, stm ] = singleDifferentialCorrection(  X_curr, designMap, constraintMap, T, tDesignFlag, simparams  )
%singleDifferentialCorrection

% Defining options for ode78
options = odeset('AbsTol',1e-12,'RelTol',1e-12);

%% Perform differential corrections
correction = 1; % Set to 1 to enable correction
stm_initial = eye(6); % STM from t_0 to t_0

% Design vector map and elements
% Initial x & z position, y velocity (and appending time T later)
designEl = find(designMap);

% Constraint vector map and elements
% Y position, x velocity, and z velocity all zero at subsequent x-z plane crossing
constraintEl = find(constraintMap);

% % Mapping of STM elements
% dfMap = constraintMap' * designMap;

while correction
    % Construct state vector with STM appended
    X = [X_curr; reshape(stm_initial, 36, 1)]; 

    % Propagate state and STM 
    [X_final, stm] = stateStmProp(X_curr, T, simparams,  stm_initial);
%     [~,X] = ode78(@(t,X) cr3bp_sFrame_nd_stm_de(X), [0,T], X, options);

    % Extract initial state vector at t=0
    X_initial = X_curr;    
    
    % Design vector - the state elements at t=0 that are being modified to 
    % try and satisfy constraints at t=T
    X_design = X_initial(designEl);
    
    % If time / Period is a design variable, append T to the end of the
    % design vector
    if tDesignFlag == 1
        X_design = [X_design; T];
        Tindex = length(X_design);
    end    
    
    % Constraint vector - the elements being minimized 
    fX_constraint = [X_final(constraintEl)];
    
    % DF matrix - portions of the STM
    % Partial of constraints wrt design variables  
    % Said another way - the sensitivity of the final constraint variables
    % to variations in the initial design variables    
    DF = stm(constraintEl,designEl);
    
    % If time is a design variable, append the time derivatives of the
    % constraint variables to the DF matrix. 
    % Incorporates the sensitivity of the constraint variables to time
    if tDesignFlag == 1
        % Calculate derivatives wrt time of state at t=T
        dXdt = cr3bp_sFrame_nd_stm_de(1, [X_final(:); stm(:)], simparams.mu);
        dXdt = dXdt(1:6);
        % Append to DF matrix
        DF = [DF, dXdt(constraintEl)];
    end
    
    % Performing corrections to minimize fX_constraint
    % X_next is the adjusted design elements (corr. to initial state, time
    % if included; more options exist)
    X_next = X_design - DF' * inv(DF*DF') * fX_constraint;

    % Re-assign X_curr and T after correction
    X_curr = zeros(6,1);
    X_curr(designEl) = X_next(1:length(designEl));
    
    if tDesignFlag == 1
        T = X_next(Tindex);
    end

    % Check if fX_constraint is sufficiently small
    disp(strcat('Norm of constraint vector after correction iteration:',num2str(norm(fX_constraint))))
    if norm(fX_constraint) < 1e-10        
        X_return = X_curr;
        correction = 0;
    end


end


end
