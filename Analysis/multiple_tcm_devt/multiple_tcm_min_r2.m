% r2 - functionalized a lot of useful parts.
% working on incorporating a third TCM.
% idea: calc next min starting from first TCM along traj. maybe include the
% first or last TCM in the calc, or both, as fixed time options.

%% For now, plot TCM history and corresponding target position dispersion
% RSS throughout the trajectory
x = simparams.x0;

[stm_i, x_i_f, x_t, stm_t, t, t_s] = createStateStmHistory(x, simparams);

% Using saved dynamics, calculate total impulsive delta V
[deltaV, deltaVs_nom] = calcDeltaV(x, x_i_f, simparams);

% Using saved dynamics, calculate the 3 sigma TCM pair
tic
[tcm_min, tcm_time, dvR3sigma_tr, dvV3sigma_tr, dvR3sigma_i, dvV3sigma_i, P_tcm1, P_tcm2, tcm_total_t, dvR3sigma_t] = tcmPair_rv(x, t, stm_t, deltaVs_nom, simparams);
toc

tcmr1_time = tcm_time;

%%  Calculate the position dispersion history 

% Calculate the target dispersion covariance and dispersion covariance
% history
% tcm_time = [];

[P_target, P_t] = calc_covariance_history(x, t, stm_t, tcm_time, simparams);






%% Plot


figure
% Plotting only the RSS of the position correction term of the TCM

plot(t, dvR3sigma_t,'DisplayName','3\sigma Position correction TCM only')
hold on
% Plotting the sum of the R and V RSS TCM
plot(t, tcm_total_t,'DisplayName','3\sigma R and V correction TCM RSS sum')

ylim([0 trimmean(tcm_total_t,5)])
xlim([0 t(end)])
legend()
title('Single TCM options along trajectory')
xlabel('TCM execution time (hrs)')
ylabel('Delta V (km/hr)')


% Plot the target position dispersion rss as a function of time
rP_t = zeros(1,length(t));
for i = 1:length(t)
    rP_t(i) = sqrt(trace(P_t(1:3,1:3,i)));
end

figure
plot(t,rP_t,'LineWidth',2)
hold on

for i = 1:length(tcm_time)
    xline(tcm_time(i))
end

title('Position dispersion RSS along trajectory with single TCM (at the min DV)')
xlim([0 t(end)])
legend('Position dispersion along trajectory','TCM time')
xlabel('Time along trajectory (hrs)')




% out of curiosity, plot the P history without any TCMs
[~, P_noTcm_t] = calc_covariance_history(x, t, stm_t, [], simparams);

rP_noTcm_t = zeros(1,length(t));
for i = 1:length(t)
    rP_noTcm_t(i) = sqrt(trace(P_noTcm_t(1:3,1:3,i)));
end

figure
plot(t, rP_noTcm_t)
title('Position dispersion RSS with no TCMs along trajectory')
xlim([0 t(end)])
xlabel('Time along trajectory (hrs)')
ylabel('Position dispersion RSS (km)')






%% Idea: show only TCM magnitudes for a single TCM that result in a target RSS that meets the constraint
% First - loop through t, calc DV and calculate P_target


tcm_dv_t = zeros(1,length(t));
tcm_dv_t(end) = 1e8;
rP_tcm_time_t = zeros(1,length(t));

tic
for i = 1:length(t) - 1
    tcm_time_i = t(i);
    [Pn_tcm_time_i, tcm_dv_t(i)] = calc_covariance_tcmdv(x, t, stm_t, tcm_time_i, simparams);
%     P_tcm_time_i = calc_covariance_history(x, t, stm_t, tcm_time_i, simparams);
    rP_tcm_time_t(i) = sqrt(trace(Pn_tcm_time_i(1:3,1:3)));
end
toc

figure
plot(t(1:end-1),rP_tcm_time_t(1:end-1));
hold on
yline(1)

xlabel('Time correction time is performed (hrs)')
ylabel('Target position dispersion rss (km)')
legend(['','1 km target position RSS constraint'])
ylim([0 max(rP_tcm_time_t)])
xlim([0 t(end)])

% Of all the final dispersion options, find the ones that meet the target
% position dispersion constraint (1 km for here)

max_target_dispersion = simparams.P_max_r; %km
rP_meets_constraint_logical = rP_tcm_time_t < max_target_dispersion;
hold on
plot(t,rP_meets_constraint_logical,'LineWidth',3)


% Find the one that is the minimum delta V
% first, find the deltaV options that meet the constraint
tcm_total_t_meets_constraint = tcm_dv_t(rP_meets_constraint_logical);
[minTcm_meets_constraint, min_meets_idx] = min(tcm_total_t_meets_constraint)
t_meets = t(rP_meets_constraint_logical);
tcm_feasibleMin_time = t_meets(min_meets_idx)

tcm_feasibleMin_idx = find(t == tcm_feasibleMin_time)

%% The above in a function - find the minimum delta V single TCM that meets the target position dispersion constraint
tic
[tcm_feasibleMin_time,tcm_feasibleMin_idx] = min_dv_tcm_meets_dispersion_constraint(x, t, stm_t, simparams);
toc

%% Second way to do the above - needs testing to see if the earliest is always the best 
% It is a bit faster

tic
[final_tcm_time,final_tcm_idx] = find_earliest_tcm_meets_Pr_constraint(x, t, stm_t ,simparams);
toc

%% Adding a third maneuver
% Currently have:
% -the cheapest initial place to perform tcmr1
% -the final time to perform tcmr3 to meet the target dispersion constraint
% Now, will search for the cheapest place to perform tcmr2 between the two

% tcmr1_time already exists
tcmr1_index = find(t==tcmr1_time);
tcmrf_time = final_tcm_time;
tcmrf_index = final_tcm_idx;
% threeTCMr_totalDV = ones(1,tcmr3_index-tcmr1_index-1)*1e8;
threeTCMr_totalDV = ones(1,length(t))*1e8;

tic
for i = 1:tcmrf_index-tcmr1_index-1
    
    threeTCMr_time = [tcmr1_time, t(i+tcmr1_index), tcmrf_time];

    [~, threeTCMr_totalDV(i+tcmr1_index)] = calc_covariance_tcmdv(x, t, stm_t, threeTCMr_time, simparams); 
    
end
toc


[minDV, minIdx] = min(threeTCMr_totalDV);
threeTCMr_DV = minDV
tcmr2_time = t(minIdx);

figure
plot(t,threeTCMr_totalDV,'LineWidth',2,'DisplayName','Total TCM \delta V')
ylim([trimmean(threeTCMr_totalDV,90)*.75 trimmean(threeTCMr_totalDV,90)*1.25])
xlim([0 t(end)])
ylabel('Total TCM Delta V (km/hr)')
xlabel('TCM 2 execution time along trajectory (hrs)')
title('Optimal 3 TCMr search')
hold on
plot(t(minIdx), threeTCMr_totalDV(minIdx),'.','MarkerSize',30,'DisplayName','Minimum total TCM \delta V')
legend()
xline([tcmr1_time, tcmrf_time])


% calculate and compare the individual magnitudes of the TCMs with 2 or 3
oneTCMr_time = tcmr1_time;
twoTCMr_time = [tcmr1_time, tcmrf_time];
threeTCMr_time = [tcmr1_time, tcmr2_time, tcmrf_time];

[~, ~, threeTCM_DVs] = calc_covariance_tcmdv(x, t, stm_t, threeTCMr_time, simparams); 

[~, twoTCMr_DV, twoTCM_DVs] = calc_covariance_tcmdv(x, t, stm_t, twoTCMr_time, simparams)

[~, ~, oneTCM_DVs] = calc_covariance_tcmdv(x, t, stm_t, oneTCMr_time, simparams); 



% Plot position covariance RSS along trajectory with two TCMrs and with three TCMrs
[~, P_2TCMr_t] = calc_covariance_history(x, t, stm_t, twoTCMr_time, simparams);
[~, P_3TCMr_t] = calc_covariance_history(x, t, stm_t, threeTCMr_time, simparams);
[~, P_1TCMr_t] = calc_covariance_history(x, t, stm_t, oneTCMr_time, simparams);

rP_1TCMr_t = zeros(1,length(t));

rP_2TCMr_t = zeros(1,length(t));
rP_3TCMr_t = zeros(1,length(t));

for i = 1:length(t)
    rP_1TCMr_t(i) = sqrt(trace(P_1TCMr_t(1:3,1:3,i)));
    rP_2TCMr_t(i) = sqrt(trace(P_2TCMr_t(1:3,1:3,i)));
    rP_3TCMr_t(i) = sqrt(trace(P_3TCMr_t(1:3,1:3,i)));
end


figure
plot(t,rP_2TCMr_t,'DisplayName','Position dispersion RSS along 2 TCMr trajectory')
hold on
plot(t,rP_3TCMr_t,'DisplayName','Position dispersion RSS along 3 TCMr trajectory')
plot(t,rP_1TCMr_t,'DisplayName','Position dispersion RSS along 1 TCMr trajectory')
xlabel('Time along trajectory (hrs)')
ylabel('Position dispersion magnitude (km)')

xline([oneTCMr_time, twoTCMr_time, threeTCMr_time])
legend()

%% Add a fourth maneuver between tcmr2 and tcmrf


tcmr2_index = find(t==tcmr2_time);
fourTCMr_totalDV = ones(1,length(t))*1e8;


for i = 1:tcmrf_index-tcmr2_index-1
    
    fourTCMr_time = [tcmr1_time, tcmr2_time, t(i+tcmr2_index), tcmrf_time];

    [~, fourTCMr_totalDV(i+tcmr2_index)] = calc_covariance_tcmdv(x, t, stm_t, fourTCMr_time, simparams); 
    
end


[minDV, minIdx] = min(fourTCMr_totalDV);

tcmr3_time = t(minIdx);
tcmr3_index = minIdx;

fourTCMr_time = [tcmr1_time, tcmr2_time, tcmr3_time, tcmrf_time];

[~, ~, fourTCM_DVs] = calc_covariance_tcmdv(x, t, stm_t, fourTCMr_time, simparams); 



figure
plot(t,fourTCMr_totalDV,'LineWidth',2,'DisplayName','Total TCM \delta V')
ylim([100 120])
xlim([0 t(end)])


ylabel('Total TCM Delta V (km/hr)')
xlabel('TCM 3 execution time along trajectory (hrs)')
title('Optimal 4 TCMr search')
hold on
plot(tcmr3_time, fourTCMr_totalDV(minIdx),'.','MarkerSize',30,'DisplayName','Minimum total TCM \delta V')
% legend()
xline(fourTCMr_time)


fourTCMr_time_best = fourTCMr_time;
fourTCMr_idx_best = [tcmr1_index, tcmr2_index, tcmr3_index, tcmrf_index];

% Get baseline delta V for comparison
[~, minDV] = calc_covariance_tcmdv(x, t, stm_t, fourTCMr_time_best, simparams); 

% The above did not result in the optimal...search for best now
tic
[fourTCMr_time_best,fourTCMr_idx_best,minDV] = tcm_index_gradient_search(x, t, stm_t, fourTCMr_time_best, fourTCMr_idx_best, 5, simparams);
toc

for i = 1:length(fourTCMr_time_best)
    plot(fourTCMr_time_best(i),minDV,'.','MarkerSize',25,'Color','Black')
end


[~, DV_gradient, fourTCM_DVs_gradient] = calc_covariance_tcmdv(x, t, stm_t, fourTCMr_time_best, simparams); 



% The brute force optimimum found
fourTCMr_idx_bruteForce = [80, 316, 419, 437];
fourTCMr_time_bruteForce = t(fourTCMr_idx_bruteForce)';
[~, DV_bruteForce, fourTCM_DVs_bruteForce] = calc_covariance_tcmdv(x, t, stm_t, fourTCMr_time_bruteForce, simparams); 



% Set baseline

fourTCMr_idx_best = [tcmr1_index, tcmr2_index, tcmr3_index, tcmrf_index];
fourTCMr_time_best = [tcmr1_time, tcmr2_time, tcmr3_time, tcmrf_time];

% Random initial selection
% fourTCMr_idx_best = sort(randi([1 length(t)],1,4));
% fourTCMr_time_best = t(fourTCMr_idx_best)';

[~, minDV] = calc_covariance_tcmdv(x, t, stm_t, fourTCMr_time_best, simparams);

startIdx = 2;
lengthMod = length(fourTCMr_time) - startIdx;
% Adding on a randomized modification method
improving = 1;
notImprovedCount = 0;
fourTCMr_idx_test = fourTCMr_idx_best;
tic
while improving
% for i = 1:100
    r1 = randi([0 1],1,lengthMod);
    r2 = randi([0 1],1,lengthMod);
    r1(logical(r2)) = r1(logical(r2)) .* -r2(logical(r2));


    fourTCMr_idx_test(startIdx:length(fourTCMr_time)-1) = fourTCMr_idx_best(startIdx:length(fourTCMr_time)-1) + r1;

    fourTCMr_time_test = t(fourTCMr_idx_test)';

    [~, testDV] = calc_covariance_tcmdv(x, t, stm_t, fourTCMr_time_test, simparams); 
    
    if testDV <= minDV
        fourTCMr_time_best = fourTCMr_time_test;
        fourTCMr_idx_best = fourTCMr_idx_test;
        minDV = testDV;
    else
        notImprovedCount = notImprovedCount + 1;
    end
    
    if notImprovedCount > 2*lengthMod^6
        improving = 0
    end


% end
end
toc
fourTCMr_DV = minDV

for i = 1:length(fourTCMr_time_best)
    plot(fourTCMr_time_best(i),minDV,'.','MarkerSize',25,'Color','Green')
end

fourTCMr_idx_test2 = [80, 317, 420, 437];
fourTCMr_time_test2 = t(fourTCMr_idx_test2)';
[~, testDV2] = calc_covariance_tcmdv(x, t, stm_t, fourTCMr_time_test2, simparams); 



%% Add a 5th maneuver

% tcmr2_index = find(t==tcmr2_time);
fiveTCMr_totalDV = ones(1,length(t))*1e8;
fiveTCMr_time = [fourTCMr_time_best(1:2), 0, fourTCMr_time_best(3:4)];
tcmr2_index_fourTCMs = fourTCMr_idx_best(2);

for i = 1:tcmrf_index-tcmr2_index_fourTCMs-1
    
%     fiveTCMr_time = [tcmr1_time, tcmr2_time, t(i+tcmr2_index_fourTCMs), tcmrf_time];
    if sum(  t(i+tcmr2_index_fourTCMs) == fiveTCMr_time  )
        % do nothing
    else
        fiveTCMr_time = sort([fourTCMr_time_best(1:2), t(i+tcmr2_index_fourTCMs), fourTCMr_time_best(3:4)]);
    
    
        [~, fiveTCMr_totalDV(i+tcmr2_index_fourTCMs)] = calc_covariance_tcmdv(x, t, stm_t, fiveTCMr_time, simparams); 
    end
    
end


[minDV, minIdx] = min(fiveTCMr_totalDV);

tcmr4_time = t(minIdx);
tcmr4_index = minIdx;

% fourTCMr_time = [tcmr1_time, tcmr2_time, tcmr3_time, tcmrf_time];
fiveTCMr_time = sort([fourTCMr_time_best(1:2), tcmr4_time, fourTCMr_time_best(3:4)]);
fiveTCMr_idx = sort( [fourTCMr_idx_best, tcmr4_index] );

[~, fiveTCM_DV, fiveTCM_DVs] = calc_covariance_tcmdv(x, t, stm_t, fiveTCMr_time, simparams); 



figure
plot(t,fiveTCMr_totalDV,'LineWidth',2,'DisplayName','Total TCM \delta V')
ylim([100 120])
xlim([0 t(end)])


ylabel('Total TCM Delta V (km/hr)')
xlabel('TCM 3 execution time along trajectory (hrs)')
title('Optimal 5 TCMr search')
hold on
plot(tcmr4_time, fiveTCMr_totalDV(minIdx),'.','MarkerSize',30,'DisplayName','Minimum total TCM \delta V')
% legend()
xline(fourTCMr_time)




[fiveTCMr_time_best,fiveTCMr_idx_best,minDV] = tcm_index_gradient_search(x, t, stm_t, fiveTCMr_time, fiveTCMr_idx, 5, simparams);

for i = 1:length(fiveTCMr_time_best)
    plot(fiveTCMr_time_best(i),minDV,'.','MarkerSize',25,'Color','Black')
end

% now the randomness improvement mod
startIdx = 2;
lengthMod = length(fiveTCMr_time) - startIdx;
% Adding on a randomized modification method
improving = 1;
notImprovedCount = 0;
fiveTCMr_idx_test = fiveTCMr_idx_best;

while improving
% for i = 1:100
    r1 = randi([0 1],1,lengthMod);
    r2 = randi([0 1],1,lengthMod);
    r1(logical(r2)) = r1(logical(r2)) .* -r2(logical(r2));


    fiveTCMr_idx_test(startIdx:length(fiveTCMr_time_best)-1) = fiveTCMr_idx_best(startIdx:length(fiveTCMr_time_best)-1) + r1;

    fiveTCMr_time_test = t(fiveTCMr_idx_test)';

    [~, testDV] = calc_covariance_tcmdv(x, t, stm_t, fiveTCMr_time_test, simparams); 
    
    if testDV <= minDV
        fiveTCMr_time_best = fiveTCMr_time_test;
        fiveTCMr_idx_best = fiveTCMr_idx_test;
        minDV = testDV;
    else
        notImprovedCount = notImprovedCount + 1;
    end
    
    if notImprovedCount > lengthMod^6
        improving = 0
    end


% end
end
fiveTCMr_DV = minDV

for i = 1:length(fiveTCMr_time_best)
    plot(fiveTCMr_time_best(i),minDV,'.','MarkerSize',25,'Color','Green')
end

%% Manual computation for verification



% 
% % Brute force test comparison -- loop over 2 parameters
% minDV = 99999;
% tcm_time_best = [0, 0, 0];
% tcmSav = zeros(length(t)-1,length(t)-1);
% 
% loopMax = cast(length(t)-1,'uint16');
% 
% tic
% for i = 1:tcmrf_index-1
%     i
%     for j = i+1:tcmrf_index-1
%         if i == j
%             % skip
%         else
% %             j
%             tcm_time_opt = [t(i), t(j), tcmrf_time];
% %             tcmSav(i,j) = two_tcmr_dv(x, t, stm_t, tcm_time_opt, simparams);
% 
% %             two_tcm_total_dv_test = two_tcmr_dv(x, t, stm_t, tcm_time_opt, simparams);
%             [~, two_tcm_total_dv_test] = calc_covariance_tcmdv(x, t, stm_t, tcm_time_opt, simparams);
%     
%             if two_tcm_total_dv_test < minDV
%                 minDV = two_tcm_total_dv_test;
%                 tcm_time_best = tcm_time_opt;
%             end
%         end
% 
%     end
% end
% toc
% 
% [~, ~, threeTCM_DVs_bruteForce] = calc_covariance_tcmdv(x, t, stm_t, tcm_time_best, simparams); 



%% Verify 4 TCMr optimal solution
% Brute force test comparison -- loop over 3 parameters
% minDV = 99999;
% tcm_time_best = [0, 0, 0, 0];


% tic
% for i = 1:tcmrf_index-1
%     i
%     for j = i+1:tcmrf_index-1
%         for k = j+1:tcmrf_index-1
% 
%             if i == j || i == k || j == k
%                 % skip
%             else
% 
%                 tcm_time_opt = [t(i), t(j), t(k), tcmrf_time];
%     
%                 [~, four_tcm_total_dv_test] = calc_covariance_tcmdv(x, t, stm_t, tcm_time_opt, simparams);
%         
%                 if four_tcm_total_dv_test < minDV
%                     minDV = four_tcm_total_dv_test;
%                     tcm_time_best = tcm_time_opt;
%                 end
%             end
%         end
% 
%     end
% end
% toc

% [~, ~, fourTCM_DVs_bruteForce] = calc_covariance_tcmdv(x, t, stm_t, tcm_time_best, simparams); 


ppp=1;

% cd('./Analysis/multiple_tcm_devt')
% % save('four_tcm_brute_force')
% load('four_tcm_brute_force')
% cd('../../')

