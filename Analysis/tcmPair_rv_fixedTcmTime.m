function [tcm_min, tcm_time, dvR3sigma_tr, dvV3sigma_tr, dvR3sigma_i, dvV3sigma_i, P_tcm1, P_tcm2, tcm_total_t] = tcmPair_rv_fixedTcmTime(x, t, stm_t, deltaVs_nom, simparams)
%tcmPair_rv Calculates the minimum TCM sum along a trajectory
% There are two TCMs required to hit a target state. The first occurs along
% the trajectory at an optimal execution time which targets a zero final
% position dispersion. The second TCM occurs at the end of the trajectory
% and corrects any final velocity dispersion.

P_initial = simparams.P_initial;
x = reshape(x,simparams.m,simparams.n);

% Extract the nominal delta V's that were passed in into separate vectors
dv1 = deltaVs_nom(:,1);
i_dv1 = dv1' / norm(dv1); % unit vector in direction of dv1
dv2 = deltaVs_nom(:,2);
i_dv2 = dv2' / norm(dv2); % unit vector in direction of dv2


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


if ~isfield(simparams,'idv_tcmV_method')
    simparams.idv_tcmV_method = 0;
end

if simparams.nom_dvctied
    % Only run through the for loop below once, at the correct index for
    % when the first maneuver occurs
    % The TCM occurs at the same time as the nominal maneuver:
    tcm_seg = simparams.maneuverSegments(simparams.maneuver_w_corr);
    tcm_time = sum(x(7,1:tcm_seg - 1));
    tcm_idx = find(tcm_time == t);

    t = tcm_time;

    stmC0 = stm_t(:,:,tcm_idx);
    stm0C = -J * stmC0' * J;
    stmNC = stmN0 * stm0C;
    T = [-inv( stmNC(1:3,4:6) ) * stmNC(1:3,1:3), -eye(3)];
    P_tcm1 = T * stmC0 * P_initial * stmC0' * T';
    


    W = [stmNC(4:6,4:6), -stmNC(4:6,1:3)];
    L = [W*T', zeros(3,3)];
    P_tcm2 = L * stmC0 * P_initial * stmC0' * L';

    dvR3sigma_tr = 3 * sqrt( trace( P_tcm1 ) );
    dvR3sigma_i = 3 * sqrt( i_dv1 * P_tcm1 * i_dv1' );

    if isfield(simparams,'idv_tcmr_method')
        simparams.idv_tcmR_method = simparams.idv_tcmr_method;
    end
    if ~isfield(simparams,'idv_tcmR_method')
        simparams.idv_tcmR_method = 0;
    end
    
    if ~simparams.idv_tcmR_method
        dvR3sigma = dvR3sigma_tr;
    else
        %%%%%%%%%%
        % compare this method with dvR3sigma and dvV3sigma
        dvR3sigma = dvR3sigma_i;
    
        %%%%%%%%%%
    end
    
    dvV3sigma_tr = 3 * sqrt( trace( P_tcm2 ) );
    dvV3sigma_i = 3 * sqrt( i_dv2 * P_tcm2 * i_dv2' );


    
    
    if ~simparams.idv_tcmV_method
        dvV3sigma = dvV3sigma_tr;
    else
        %%%%%%%%%%
        % compare this method with dvR3sigma and dvV3sigma
        dvV3sigma = dvV3sigma_i;
    
        %%%%%%%%%%
    end




    tcm_total_t = dvR3sigma + dvV3sigma;
    tcm_min = tcm_total_t;
    dvR3sigma_return = dvR3sigma;

    dvV3sigma_return = dvV3sigma;

else
    % Otherwise, run through the entire trajectory to find the lowest TCM

    % Allocation as needed
    dvR3sigma = zeros(size(t,1),1);
    dvV3sigma = zeros(size(t,1),1);
    dvV3sigma_tr = zeros(size(t,1),1);
    dvV3sigma_i = zeros(size(t,1),1);
    
    
    % Calculate the minimum TCM pair
    for i = 1:size(t,1)
    
    %     t_c = t(i);
    %     x_c = x_t(i,:)';
        stmC0 = stm_t(:,:,i); %%%%%%%% BUG
        stm0C = -J * stmC0' * J; % symplectic inverse 
        stmNC = stmN0 * stm0C;
    
        T = [-inv( stmNC(1:3,4:6) ) * stmNC(1:3,1:3), -eye(3)];
        P_tcm1 = T * stmC0 * P_initial * stmC0' * T';

        
    
        % Magnitude of correction to final position dispersion
        dvR3sigma(i) = 3 * sqrt( trace( P_tcm1 ) );

        N = [zeros(3,6); T];
        IN = eye(6) + N;
    
        % Magnitude of correction to final velocity dispersion
        W = [stmNC(4:6,4:6), -stmNC(4:6,1:3)];
        L = [W*T', zeros(3,3)];
%         P_tcm2 = L * IN * stmC0 * P_initial * stmC0' * IN' * L';
        P_tcm2 = L  * stmC0 * P_initial * stmC0' * L';

        %%%%%%%% test/compare the following 
        
        dvV3sigma_tr(i) = 3 * sqrt( trace( P_tcm2 ) ); %%%% old way, assumes gaussian
%         dvV3sigma_i(i) = 3 * sqrt( i_dv2 * P_tcm2 * i_dv2' ); %%%% new way, needs testing...is underestimating a bit in some circumnstances

    
    end

    if ~simparams.idv_tcmV_method
        dvV3sigma = dvV3sigma_tr;
    else
        dvV3sigma = dvV3sigma_i;
    end

    
    tcm_total_t = dvR3sigma + dvV3sigma;
    [tcm_min, tcm_min_idx] = min(tcm_total_t);
    
    tcm_time = t(tcm_min_idx);

    dvR3sigma_tr = dvR3sigma(tcm_min_idx);
    dvR3sigma_i = nan;
    dvV3sigma_tr = dvV3sigma_tr(tcm_min_idx);
    dvV3sigma_i = dvV3sigma_i(tcm_min_idx);
    


end





end