

%% For now, plot TCM history and corresponding target position dispersion
% RSS throughout the trajectory
x = simparams.x0;

[stm_i, x_i_f, x_t, stm_t, t, t_s] = createStateStmHistory(x, simparams);

% Using saved dynamics, calculate total impulsive delta V
[deltaV, deltaVs_nom] = calcDeltaV(x, x_i_f, simparams);

% Using saved dynamics, calculate the 3 sigma TCM pair

[tcm_min, tcm_time, dvR3sigma_tr, dvV3sigma_tr, dvR3sigma_i, dvV3sigma_i, P_tcm1, P_tcm2, tcm_total_t, dvR3sigma_t] = tcmPair_rv(x, t, stm_t, deltaVs_nom, simparams);

f_idx = round(length(t)*.9);



%%  Calculate the position dispersion history 

% Calculate the target dispersion covariance
% tcm_time = [];

[P_target, P_t] = calc_covariance_history(x, t, stm_t, tcm_time, simparams);


% Calculate the target dispersion covariance history




%% Plot

% Plotting ONLY the position correction term
figure
plot(t(1:f_idx), dvR3sigma_t(1:f_idx),'DisplayName','3\sigma Position correction TCM only')
hold on
% Plotting the sum of the the R and V 3 sigma corrections
plot(t(1:f_idx), tcm_total_t(1:f_idx),'DisplayName','3\sigma R and V correction TCM RSS sum')

ylim([0 1000])
legend()


% Plotting ONLY the position correction term
figure
plot(t, dvR3sigma_t,'DisplayName','3\sigma Position correction TCM only')
hold on
% Plotting the sum of the the R and V 3 sigma corrections
plot(t, tcm_total_t,'DisplayName','3\sigma R and V correction TCM RSS sum')

ylim([0 1000])
xlim([0 t(end)])
legend()


% Plot the position dispersion rss as a function of time
rP_t = zeros(1,length(t));
for i = 1:length(t)
    rP_t(i) = sqrt(trace(P_t(1:3,1:3,i)));
end

figure
plot(t,rP_t)
hold on
for i = 1:length(tcm_time)
    xline(tcm_time(i))
end

title('Position dispersion RSS with single TCM (min)')


% out of curiosity, plot the P history without any TCMs
[~, P_noTcm_t] = calc_covariance_history(x, t, stm_t, [], simparams);

rP_noTcm_t = zeros(1,length(t));
for i = 1:length(t)
    rP_noTcm_t(i) = sqrt(trace(P_noTcm_t(1:3,1:3,i)));
end

figure
plot(t, rP_noTcm_t)
title('Position dispersion RSS with no TCMs')






%% Idea: show only TCM magnitudes for a single TCM that result in a target RSS that meets the constraint
% First - loop through t, calc DV (tcm1 and tcmv (already exists in tcm_total_t)) and calculate P_target

tcm_dv_t = zeros(1,length(t));
rP_tcm_time_t = zeros(1,length(t));

for i = 1:length(t) - 1
    tcm_time_i = t(i);
    P_tcm_time_i = calc_covariance_history(x, t, stm_t, tcm_time_i, simparams);
    rP_tcm_time_t(i) = sqrt(trace(P_tcm_time_i(1:3,1:3)));


end

figure
plot(t(1:end-1),rP_tcm_time_t(1:end-1));
hold on
yline(1)

xlabel('correction time (hrs)')
ylabel('target position dispersion rss (km)')
legend(['','1 km target position RSS constraint'])

% Of all the final dispersion options, find the ones that meet the target
% position dispersion constraint (1 km for here)

max_target_dispersion = 1; %km
rP_meets_constraint_logical = rP_tcm_time_t < max_target_dispersion;
hold on
plot(t,rP_meets_constraint_logical,'LineWidth',3)

% Find the one that is the minimum delta V
% first, find the deltaV options that meet the constraint
tcm_total_t_meets_constraint = tcm_total_t(rP_meets_constraint_logical);
[minTcm_meets_constraint, min_meets_idx] = min(tcm_total_t_meets_constraint)
t_meets = t(rP_meets_constraint_logical);
tcm_feasibleMin_time = t_meets(min_meets_idx)

tcm_feasibleMin_idx = find(t == tcm_feasibleMin_time)

%% Now find the lowest DV place to incorporate a second TCM
% the option set is all the times earlier than tcm_feasibleMin_time

tcm2_time_options = t( t < tcm_feasibleMin_time );

m = simparams.m;
n = simparams.n;
x = reshape(x,m,n);

% calc stmN0
if simparams.target_final_maneuver
    final_seg = simparams.maneuverSegments(end);
    target_time = sum(x(7,1:final_seg - 1));
    target_idx = find(target_time == t) ;
    stmN0 = stm_t(:,:,target_idx);

else
    stmN0 = stm_t(:,:,end);
end



tic
for i = 1:length(tcm2_time_options)
    tcm_time_opt = [tcm2_time_options(i), tcm_feasibleMin_time];

    two_tcm_total_dv(i) = two_tcmr_dv(x, t, stm_t, tcm_time_opt, simparams);


    



end
toc



figure
plot(tcm2_time_options,two_tcm_total_dv)
ylim([0 1000])




[tcm1_total_min_dv, tcm1_time_idx] = min(two_tcm_total_dv);
tcm1_min_time = t(t == tcm2_time_options(tcm1_time_idx));




% Brute force test comparison
minDV = 99999;
tcm_time_best = [0, 0];
tcmSav = zeros(length(t)-1,length(t)-1);

loopMax = cast(length(t)-1,'uint16');

tic
for i = 1:loopMax
    i
    for j = tcm_feasibleMin_idx:loopMax
        if i == j
            % skip
        else
%             j
            tcm_time_opt = [t(i), t(j)];
%             tcmSav(i,j) = two_tcmr_dv(x, t, stm_t, tcm_time_opt, simparams);

%             two_tcm_total_dv_test = two_tcmr_dv(x, t, stm_t, tcm_time_opt, simparams);
            [~, two_tcm_total_dv_test] = calc_covariance_tcmdv(x, t, stm_t, tcm_time_opt, simparams);
    
            if two_tcm_total_dv_test < minDV
                minDV = two_tcm_total_dv_test;
                tcm_time_best = tcm_time_opt;
            end
        end

    end
end
toc


ppp=1;
