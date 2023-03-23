dx = 1e-5;

x = simparams.x0;
x=x_opt;

% Determining the partial derivative of tcm_time wrt problem params

% First step - unperturbed propagation
[stm_i, x_i_f, x_t, stm_t, t] = createStateStmHistory(x, simparams);
[deltaV, deltaVs_nom] = calcDeltaV(x, x_i_f, simparams);
[tcm_3sigma,tcm_time, dvR3sigma_tr, dvV3sigma_tr] = tcmPair_rv(x, t, stm_t, deltaVs_nom, simparams);

% Empty structure to store dTcm_time_dx
dTcm_time_dx = zeros(1,length(x));

% Empty structure to store dJ_fixTcm_dX
dJ_fixTcm_dX = zeros(1,length(x));

for j = 1:length(x(:)) 
    xdx = x;
    xdx(j) = x(j) + dx;

    % Propagate perturbed traj (_p)

    [stm_i_p, x_i_f_p, x_t_p, stm_t_p, t_p] = createStateStmHistory(xdx, simparams);
    [deltaV_p, deltaVs_nom_p] = calcDeltaV(xdx, x_i_f_p, simparams);
    [tcm_3sigma_p,tcm_time_p, dvR3sigma_tr_p, dvV3sigma_tr_p] = tcmPair_rv(xdx, t_p, stm_t_p, deltaVs_nom_p, simparams);

    % Calculate and store sensitivity of tcm time wrt to xdx
    dTcm_time_dx(j) = (tcm_time_p - tcm_time)/dx;
        %%%% There is some/a tiny bit...mostly zeros though.
        % Largest number appearing in some time sensitivity terms (close to
        % 1)...is there a 1 needed somewhere in some logic case in the
        % analytical gradient?
    




    % Perform another test - verify analytical gradient solution without
    % the optimal tcm time step. See if the small discrepancies are gone.


    




    
end
dTcm_time_dx'



pp=1;

