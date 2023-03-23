clear;
clear global;
close all; 
clc;
format longg;
cd('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust_2body_tcmerror');
addpath(genpath('./'));



init_fn = './init_traj_files/init_simparams_cr3bp_leo_to_llo';

% init_fn = './init_traj_files/initialize_simulation_parameters_leo28_to_geo0';



run(init_fn);

mu = simparams.mu;

% simparams = generateInitialXfer(simparams);
% ndTime2days = 1;
% ndDist2km = 1;
% ndVel2kms = 1;


x = simparams.x0;
x = reshape(x,simparams.m,simparams.n);


[~, ~, x_t, stm_t, t, t_s] = createStateStmHistory(simparams.x0, simparams);
[tcm_opt_times,tcm_idx,tcm_opt_dv] = opt_multiple_tcm(simparams.x0, t, stm_t, simparams);


if simparams.target_final_maneuver
    final_seg = simparams.maneuverSegments(end);
    target_time = sum(x(7,1:final_seg - 1));
    target_idx = find(target_time == t) ;
    stmN0 = stm_t(:,:,target_idx);
    % Abbreviate t and stm_t to go only until the target (not the end of
    % the trajectory)
    t = t(1:target_idx);
    stm_t = stm_t(:,:,1:target_idx);
else
    stmN0 = stm_t(:,:,end);    
end



% Plot the nominal trajectory
figure;
% plotMultiSegTraj(simparams.x0, x_t0, t_s0, simparams, tcm_idx0);
plotMultiSegTraj(simparams.x0, x_t, t_s, simparams);
axis equal




%%
%min_dv_tcm_meets_dispersion_constraint Finds all feasible TCM options that
%meet the target dispersion constraint, then selects the lowest DV option
%of those and returns the time and time index.

tcm_dv_t = zeros(1,length(t));
tcm_dv_t(end) = 1e8; % the loop below skips the final time index...a correction once at the final time doesn't make sense.
rP_tcm_time_t = zeros(1,length(t));


% Loop through the time options, save the target position dispersion and
% the TCM delta V RSS if execute at that time
for i = 1:length(t) - 1
    tcm_time_i = t(i);
    [Pn_tcm_time_i, tcm_dv_t(i)] = calc_covariance_tcmdv(x, t, stm_t, tcm_time_i, simparams);
    rP_tcm_time_t(i) = sqrt(trace(Pn_tcm_time_i(1:3,1:3)));
end



%%%%%% testing 

% [Pn_test, tcm_dv_test] = calc_covariance_tcmdv(x, t, stm_t, t(764), simparams);


%%%%%%%%




% Of all the final dispersion options, find the ones that meet the target
% position dispersion constraint (1 km for here)

rP_meets_constraint_logical = rP_tcm_time_t <= simparams.P_max_r;

rP_meets_constraint_logical(end) = false;

% Find the one that is the minimum delta V
% first, find the deltaV options that meet the constraint
tcm_total_t_meets_constraint = tcm_dv_t(rP_meets_constraint_logical);
% the find the minimum, and the index of the minimum
[minTcmDV_meets_constraint, min_meets_idx] = min(tcm_total_t_meets_constraint);
% find the corresponding time index of all tcm options that meet the constraint
t_meets = t(rP_meets_constraint_logical);
% get the time from the feasible set index and feasible set time array
tcm_feasibleMin_time = t_meets(min_meets_idx);
tcm_feasibleMin_idx = find(t == tcm_feasibleMin_time);


%%




%%%%%% DEBUGGING:

P_t = tmult(stm_t, tmult(simparams.P_initial, stm_t, [0 1]));

for i = 1:size(P_t,3)
    Pmag_t(i) = sqrt(trace(P_t(1:3,1:3,i)));
end
figure;
plot(t*ndTime2days,Pmag_t*ndDist2km,'LineWidth',2)
xlim([t(1)*ndTime2days t(end)*ndTime2days]);
xlabel('Time along nominal trajectory (days)')
ylabel('Position dispersion RSS (km)')
title('Uncorrected position dispersion')


% Calculate dmagPdt: rate of change of P over time
for i = 2:length(t)
    dmagP_t(i) = (Pmag_t(i) - Pmag_t(i-1))/(t(i) - t(i-1));

end

figure
plot(t*ndTime2days, dmagP_t,'LineWidth',2)
xlim([t(1)*ndTime2days t(end)*ndTime2days]);
title('Rate of change of position dispersion RSS')
xlabel('Time along nominal trajectory (days)')
ylabel('Position dispersion RSS rate (km/s)')


figure
plot(t*ndTime2days,rP_tcm_time_t*ndDist2km,'LineWidth',2,'DisplayName','Pr_n')
xlim([t(1)*ndTime2days t(end)*ndTime2days]);
hold on
yline(simparams.P_max_r * ndDist2km,'DisplayName','Max target Pr constraint')
plot(t*ndTime2days,rP_meets_constraint_logical * simparams.P_max_r * ndDist2km,'LineWidth',2,'DisplayName','Meets Pr_n constraint flag')
title('Target position dispersion as a function of TCM execution time (w/error)')
xlabel('Time along nominal trajectory (days)')
ylabel('Position dispersion RSS (km)')
legend()




tcm_dv_t(end) = 5;
figure
plot(t*ndTime2days,tcm_dv_t*ndVel2kms,'LineWidth',2,'DisplayName','TCM \deltaV')
hold on
plot(t*ndTime2days,rP_meets_constraint_logical * .1,'LineWidth',2,'DisplayName','Meets Pr_n constraint flag')

plot(tcm_feasibleMin_time * ndTime2days, minTcmDV_meets_constraint*ndVel2kms,'.','MarkerSize',25,'DisplayName','Min feasible TCM \deltaV')


xlim([t(1)*ndTime2days t(end)*ndTime2days]);
ylim([0 .2])
xlabel('Time along nominal trajectory (days)')
ylabel('TCM RSS (km/s)')
legend()








assert(0)



%% Options for adding a second TCM
for i = 1:length(t) - 1
    tcm_time_2_i = t(i);
    if tcm_time_2_i ~= tcm_feasibleMin_time
        [Pn_tcm_time_i, tcm_2_dv_t(i)] = calc_covariance_tcmdv(x, t, stm_t, sort([tcm_feasibleMin_time, tcm_time_2_i]), simparams);
        rP_2_tcm_time_t(i) = sqrt(trace(Pn_tcm_time_i(1:3,1:3)));
    else
        tcm_2_dv_t(i) = tcm_dv_t(i);
        rP_2_tcm_time_t(i) = rP_tcm_time_t(i);
    end


end


figure
plot(t(1:end-1)*ndTime2days,rP_2_tcm_time_t*ndDist2km,'LineWidth',2,'DisplayName','Pr_n')
xlim([t(1)*ndTime2days t(end)*ndTime2days]);
hold on
yline(simparams.P_max_r * ndDist2km,'DisplayName','Max target Pr constraint')
% plot(t*ndTime2days,rP_meets_constraint_logical * simparams.P_max_r * ndDist2km,'LineWidth',2,'DisplayName','Meets Pr_n constraint flag')
plot(tcm_feasibleMin_time * ndTime2days, rP_2_tcm_time_t(tcm_feasibleMin_idx)*ndDist2km,'.','MarkerSize',25,'DisplayName','Fixed TCM execution time')
title('Target position dispersion as a function of 2nd TCM execution time (w/error)')
xlabel('Time along nominal trajectory the 2nd TCM is applied (days)')
ylabel('Target Position dispersion RSS (km)')
legend()






tcm_2_dv_t(end) = 5;
figure
plot(t(1:end-1)*ndTime2days,tcm_2_dv_t*ndVel2kms,'LineWidth',2,'DisplayName','TCM \deltaV')
hold on
% plot(t*ndTime2days,rP_meets_constraint_logical * .1,'LineWidth',2,'DisplayName','Meets Pr_n constraint flag')

% plot(tcm_feasibleMin_time * ndTime2days, minTcmDV_meets_constraint*ndVel2kms,'.','MarkerSize',25,'DisplayName','Min feasible TCM \deltaV')






[minDv2, minDv2_idx] = min(tcm_2_dv_t);

plot(t(minDv2_idx) * ndTime2days, minDv2*ndVel2kms,'.','MarkerSize',25,'DisplayName','Min direct search for 2nd TCM execution time')


xlim([t(1)*ndTime2days t(end)*ndTime2days]);
ylim([0 .2])
xlabel('Time along nominal trajectory (days)')
ylabel('TCM RSS (km/s)')
legend()



% Pick the minimum of the options that meet the target position dispersion
% constraint
rP_2_meets_constraint_logical = rP_2_tcm_time_t <= simparams.P_max_r;

rP_2_meets_constraint_logical(end) = false;



tcm_2_total_t_meets_constraint = tcm_2_dv_t(rP_2_meets_constraint_logical);
% the find the minimum, and the index of the minimum
[minTcm2DV_meets_constraint, min2_meets_idx] = min(tcm_2_total_t_meets_constraint);
% find the corresponding time index of all tcm options that meet the constraint
t2_meets = t(rP_2_meets_constraint_logical);
% get the time from the feasible set index and feasible set time array
tcm2_feasibleMin_time = t2_meets(min2_meets_idx);
tcm2_feasibleMin_idx = find(t == tcm2_feasibleMin_time);


% 2 TCM guess from sequential direct search
tcm_time_2dv = sort([tcm2_feasibleMin_time, tcm_feasibleMin_time]);
tcm_idx_2dv = sort([tcm2_feasibleMin_idx, tcm_feasibleMin_idx]);


% Test if gradient might work
tcm_time_2dv_test = sort([tcm2_feasibleMin_time, t(tcm_feasibleMin_idx-1)]);

[~, test_dv] = calc_covariance_tcmdv(x, t, stm_t, tcm_time_2dv_test, simparams);
minTcm2DV_meets_constraint - test_dv



% Perform gradient search method
[tcm_idx_test, tcm_time_test, minDV] = tcm_index_gradient_vector_search(x, t, stm_t, tcm_idx_2dv, 1, simparams);
% Didn't go all the way to the min (the first TCM shouldn't move but the
% second TCM didn't move early enough to find the min)

% Look at options by holding TCM 1 fixed
for i = 1:length(t) - 1
    tcm_time_2_i = t(i);
    if tcm_time_2_i ~= t(1)
        [Pn_tcm_time_i, tcm_m2_dv_t(i)] = calc_covariance_tcmdv(x, t, stm_t, sort([t(1), tcm_time_2_i]), simparams);
        rP_m2_tcm_time_t(i) = sqrt(trace(Pn_tcm_time_i(1:3,1:3)));
    else
        tcm_m2_dv_t(i) = tcm_dv_t(i);
        rP_m2_tcm_time_t(i) = rP_tcm_time_t(i);
    end


end








figure
plot(t(1:end-1)*ndTime2days,rP_m2_tcm_time_t*ndDist2km,'LineWidth',2,'DisplayName','Pr_n')
xlim([t(1)*ndTime2days t(end)*ndTime2days]);
hold on
yline(simparams.P_max_r * ndDist2km,'DisplayName','Max target Pr constraint')
% plot(t*ndTime2days,rP_meets_constraint_logical * simparams.P_max_r * ndDist2km,'LineWidth',2,'DisplayName','Meets Pr_n constraint flag')
% plot(tcm_feasibleMin_time * ndTime2days, rP_m2_tcm_time_t(tcm_feasibleMin_idx)*ndDist2km,'.','MarkerSize',25,'DisplayName','Fixed TCM execution time')
title('Target position dispersion as a function of 2nd TCM execution time (w/error)')
xlabel('Time along nominal trajectory the 2nd TCM is applied (days)')
ylabel('Target Position dispersion RSS (km)')
legend()






tcm_m2_dv_t(end) = 5;
figure
plot(t(1:end-1)*ndTime2days,tcm_m2_dv_t*ndVel2kms,'LineWidth',2,'DisplayName','TCM \deltaV')
hold on
% plot(t*ndTime2days,rP_meets_constraint_logical * .1,'LineWidth',2,'DisplayName','Meets Pr_n constraint flag')

% plot(tcm_feasibleMin_time * ndTime2days, minTcmDV_meets_constraint*ndVel2kms,'.','MarkerSize',25,'DisplayName','Min feasible TCM \deltaV')






[minDvm2, minDvm2_idx] = min(tcm_m2_dv_t);

plot(t(minDvm2_idx) * ndTime2days, minDvm2*ndVel2kms,'.','MarkerSize',25,'DisplayName','Min direct search for 2nd TCM execution time')


xlim([t(1)*ndTime2days t(end)*ndTime2days]);
ylim([0 .2])
xlabel('Time along nominal trajectory (days)')
ylabel('TCM RSS (km/s)')
legend()




%%
% state position and velocity components
% velocity
figure
plot(t,x_t(1:target_idx,4))
hold on
plot(t,x_t(1:target_idx,5))
plot(t,x_t(1:target_idx,6))



%%




tcm_time_2_best = [1 2];
tcm_time_2_best_mag = 1e8;

% brute force calculate 2 TCM solution that meets the covariance constraint
for i = 1:length(t)
    for j = 1:length(t)
        if t(i) ~= t(j)
            [Pn_test, tcm_dv_test] = calc_covariance_tcmdv(x, t, stm_t, sort([t(i) t(j)]), simparams);
            rP_test_mag = sqrt(trace(Pn_test(1:3,1:3)));

            if tcm_dv_test < tcm_time_2_best_mag && rP_test_mag <= simparams.P_max_r
                tcm_time_2_best = sort([t(i) t(j)]);
                tcm_time_2_best_mag = tcm_dv_test;
            end

        end

    end
end

tcm_time_2_best = [1 2];
tcm_time_2_best_mag = 1e8;


% brute force 3


%% plot dispersion ellipsoids

%%%%% COME BACK HERE AND WORK ON ELLIPSOIDS!
