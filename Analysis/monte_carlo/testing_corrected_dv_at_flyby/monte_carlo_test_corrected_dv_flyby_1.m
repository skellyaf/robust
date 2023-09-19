clear;
clear global;
close all; 
clc;
format longg;
addpath(genpath('../../../'));

savename = ['testing_dv_stats_nri'];
saveOutput = true; % bool for saving the output or not, true or false
saveVideo = false;


%% Create simulation parameters structure by running initialization script

% init_fn = './init_traj_files/init_simparams_cr3bp_leo_lloflyby_nri_3dv';
init_fn = '../../../init_traj_files/init_simparams_cr3bp_leoinclined_lloflyby_nri_3dv';

% init_fn = './init_traj_files/init_simparams_cr3bp_leoinclined_lloflyby_dro_3dv';
% init_fn = './init_traj_files/init_simparams_cr3bp_leoinclined_lloflyby_dro2_3dv';
% init_fn = './init_traj_files/init_simparams_cr3bp_leoinclined_lloflyby_dro3_3dv';
% init_fn = './init_traj_files/init_simparams_cr3bp_leoinclined_lloflyby_dro4_3dv';

% init_fn = './init_traj_files/init_simparams_cr3bp_leoinclined_lloflyby_dro3_x0detOpt_3dv';
options = optimoptions(@fsolve)
options.StepTolerance = eps;
options.OptimalityTolerance = 1e-20;
options.FiniteDifferenceStepSize = eps*1e3;
options.FunctionTolerance = 1e-20;
options.Display = 'none';
options.MaxFunctionEvaluations = 1000;



% init_fn = './init_traj_files/init_simparams_cr3bp_bigHeo_to_mlo_3dv';
% init_fn = './init_traj_files/init_simparams_cr3bp_flex0_leo_lloflyby_nri_3dv';



run(init_fn);

%% ONLY If using a previous reference trajectory as the initial guess:
% load('.\init_traj_files\initial_guesses\tli_llo_deterministic_opt.mat')
% load('.\init_traj_files\initial_guesses\polar_llo_to_nrho_apolune.mat')

% load('eed_leo_planar_13day.mat');
load('nri_det_opt.mat');
% load('nri_planar_det_opt.mat');
% load('leo_plf_dro3_detOpt.mat')
simparams.x0 = x_opt;


%% Initialize transfer segments
% simparams = generateInitialXfer(simparams); % only for 2bp

%% View problem setup

% [~, x_i_f0, x_t0, stm_t0, t0, t_s0] = createStateStmHistory(simparams.x0, simparams);
[stm_i0, stt_i0, x_i_f0, x_t0, stm_t0, stt_t_i, t0, t_s0, stm_t_i0]  = createStateStmSttHistory(simparams.x0, simparams);

% Calculate total impulsive delta V for initial guess trajectory
[deltaV0, deltaVs_nom0] = calcDeltaV(simparams.x0,x_i_f0,simparams);
% [tcm_min, tcm_time, tcm_r, tcm_v, ~, ~, ~, ~, tcm_total_t] = tcmPair_rv(x_opt, t, stm_t, deltaVs_nom, simparams);

tic
[tcm_time0, tcm_idx0, min_tcm_dv0, P_i_minus, P_i_plus, tcm_dv_each0] = opt_multiple_tcm(simparams.x0, deltaVs_nom0, t0, t_s0, stm_t0, stm_t_i0, simparams); % inputs: x, t, t_s, stm_t, stm_t_i, simparams
toc

totalDV0 = deltaV0 + 3*min_tcm_dv0

figure;
if simparams.perform_correction    
%     [~,tcm_idx0] = opt_multiple_tcm(simparams.x0, deltaVs_nom0, t0, t_s0, stm_t0, stm_t_i0, simparams);
    plotMultiSegTraj(simparams.x0, x_t0, t_s0, simparams, tcm_idx0);
else
    plotMultiSegTraj(simparams.x0, x_t0, t_s0, simparams);
end
axis equal;

x = reshape(simparams.x0,simparams.m,simparams.n);

%% Perform monte carlo
nRuns = 500;

% Data storage
dv1_correction_i = zeros(3,nRuns);
dv2_correction_i = zeros(3,nRuns);
dv2_corr_norm_diff_i = zeros(1,nRuns);
tcm1_corrected_i = zeros(3,nRuns);
tcm2_corrected_i = zeros(3,nRuns);
tcm3_corrected_i = zeros(3,nRuns);

dv2_position_i = zeros(3,nRuns);



% Precompute - things that won't change each sim
x0 = x(1:6)';
% Get events
[event_times, event_indicator] = define_events_v2(simparams.x0, t0, tcm_time0, simparams);
event_idx_logical = logical(sum(t0'==event_times', 1));    
event_idxs = find(event_idx_logical);

dv_event_idx = find(event_indicator==3,2);
dv1_event_idx = dv_event_idx(1);
assert(dv1_event_idx == 1); % placeholder in case I forget to change this in a different simulation that doesn't abide
dv1_time = event_times(dv1_event_idx);
dv1_idx = event_idxs(1);



% Add the first delta V, corrected, with error
dv2_event_idx = dv_event_idx(2);
dv2_time = event_times(dv2_event_idx);
dv2_idx = event_idxs(dv2_event_idx);
dv2_seg = simparams.maneuverSegments(2);
dv3_seg = simparams.maneuverSegments(3);

dv2_dv1_time = event_times(dv2_event_idx) - event_times(dv1_event_idx);
nom_dv_1 = x(4:6,2) - x_i_f0(4:6,1);

nom_dv_2 = x(4:6,dv2_seg) - x_i_f0(4:6,dv2_seg - 1);


dv3_time = sum(x(7,1:dv3_seg-1));



t_tcm_1 = event_times(2);
tcm_1_idx = event_idxs(2);  

t_tcm_2 = event_times(3);
tcm_2_idx = event_idxs(3);
t_tcm_3 = event_times(4);
tcm_3_idx = event_idxs(4);

stm_dv2_dv1 = dynCellCombine(t0, t_s0, dv1_idx, dv2_idx, simparams, stm_t_i0);
T_dv2_dv1 = [-inv( stm_dv2_dv1(1:3,4:6) ) * stm_dv2_dv1(1:3,1:3), -eye(3)];


stm_dv2_tcm1 = dynCellCombine(t0, t_s0, tcm_1_idx, dv2_idx, simparams, stm_t_i0);
T1 = [-inv( stm_dv2_tcm1(1:3,4:6) ) * stm_dv2_tcm1(1:3,1:3), -eye(3)];

stm_dv2_tcm2 = dynCellCombine(t0, t_s0, tcm_2_idx, dv2_idx, simparams, stm_t_i0);
T2 = [-inv( stm_dv2_tcm2(1:3,4:6) ) * stm_dv2_tcm2(1:3,1:3), -eye(3)];

stm_dv2_tcm3 = dynCellCombine(t0, t_s0, tcm_3_idx, dv2_idx, simparams, stm_t_i0);
T3 = [-inv( stm_dv2_tcm3(1:3,4:6) ) * stm_dv2_tcm3(1:3,1:3), -eye(3)];


stm_dv3_dv2 = dynCellCombine(t0, t_s0, tcm_3_idx, dv2_idx, simparams, stm_t_i0);
T_dv3_dv2 = [-inv( stm_dv3_dv2(1:3,4:6) ) * stm_dv3_dv2(1:3,1:3), -eye(3)];

badidx = false(1,nRuns);
saveexitflag1 = zeros(1,nRuns);
saveexitflag2 = zeros(1,nRuns);
saveexitflag3 = zeros(1,nRuns);
saveexitflag4 = zeros(1,nRuns);
saveexitflag5 = zeros(1,nRuns);
x0_save = zeros(6,nRuns);

parfor i = 1:nRuns
% for i = 1:nRuns
    i
    % Synthesize initial error state

    % No rotation required since it is aligned with the coordinate frame to start
    x0_r = x0 + randn(6,1).*diag(simparams.P_initial); %_r indicates it is a single random state
    x0_save(:,i) = x0_r;

    

    % Propagate to first nominal maneuver
    x_dv1_minus_r = stateProp(x0_r, dv1_time, simparams); % the random state prior to the first nominal maneuver

    

    % Compute corrected dv1
    [nom_dv_1_corrected, ~, exitflag, output] = fsolve(  @(nom_dv)fsolve_tcm_target_function2(nom_dv, x_dv1_minus_r, x(1:6,dv2_seg), dv2_dv1_time, simparams), nom_dv_1, options  );
    saveexitflag1(i) = exitflag;
    if exitflag ~=3
        badidx(i) = true;
        
        ppp=1;
    end
    dv1_correction_i(:,i) = nom_dv_1_corrected - nom_dv_1;
    dv1_corr_norm_diff_i(i) = norm(nom_dv_1_corrected) - norm(nom_dv_1);

    % Add DV1 and error
    x_dv1_plus_r = x_dv1_minus_r + [zeros(3,1); nom_dv_1_corrected + randn(3,1) .* simparams.sig_dv_error];

    % Propagate to the first TCM
    x_tcm_1_minus = stateProp(x_dv1_plus_r, t_tcm_1 - event_times(dv1_event_idx), simparams);
    
    % Estimate TCM from dispersion
    dx_tcm_1 = x_tcm_1_minus - x_t0(t0==t_tcm_1,:)';
    tcm1_dv = T1 * dx_tcm_1;    

    % Calculate TCM 1
    [tcm_1_corrected, ~, exitflag, output] = fsolve(  @(tcm_dv)fsolve_tcm_target_function2(tcm_dv, x_tcm_1_minus, x(1:6,dv2_seg), event_times(dv2_event_idx) - t_tcm_1, simparams), tcm1_dv, options  );
    saveexitflag2(i) = exitflag;
    if exitflag ~=3
        badidx(i) = true;
        
        ppp=1;
    end

    tcm_1_corrected_i(:,i) = tcm_1_corrected;
    x_tcm_1_plus_r = x_tcm_1_minus + [zeros(3,1); tcm_1_corrected + randn(3,1) .* simparams.sig_tcm_error];




    % Propagate to TCM 2
    x_tcm_2_minus = stateProp(x_tcm_1_plus_r, t_tcm_2 - t_tcm_1, simparams);
    
    % Estimate TCM 2 from dispersion
    dx_tcm_2 = x_tcm_2_minus - x_t0(t0==t_tcm_2,:)';
    tcm2_dv = T2 * dx_tcm_2;    

    % Calculate TCM 2
    [tcm_2_corrected, ~, exitflag, output] = fsolve(  @(tcm_dv)fsolve_tcm_target_function2(tcm_dv, x_tcm_2_minus, x(1:6,dv2_seg), event_times(dv2_event_idx) - t_tcm_2, simparams), tcm2_dv, options  );
    saveexitflag3(i) = exitflag;
    if exitflag ~=3
        badidx(i) = true;
        
        ppp=1;
    end
    tcm_2_corrected_i(:,i) = tcm_2_corrected;
    x_tcm_2_plus_r = x_tcm_2_minus + [zeros(3,1); tcm_2_corrected + randn(3,1) .* simparams.sig_tcm_error];





    % Propagate to TCM 3
    x_tcm_3_minus = stateProp(x_tcm_2_plus_r, t_tcm_3 - t_tcm_2, simparams);
    
    % Estimate TCM 3 from dispersion
    dx_tcm_3 = x_tcm_3_minus - x_t0(t0==t_tcm_3,:)';
    tcm3_dv = T3 * dx_tcm_3;    

    % Calculate TCM 3
    [tcm_3_corrected, ~, exitflag, output] = fsolve(  @(tcm_dv)fsolve_tcm_target_function2(tcm_dv, x_tcm_3_minus, x(1:6,dv2_seg), event_times(dv2_event_idx) - t_tcm_3, simparams), tcm3_dv, options  );
    saveexitflag4(i) = exitflag;
    if exitflag ~=3
        badidx(i) = true;
        
        ppp=1;
    end
    tcm_3_corrected_i(:,i) = tcm_3_corrected;
    x_tcm_3_plus_r = x_tcm_3_minus + [zeros(3,1); tcm_3_corrected + randn(3,1) .* simparams.sig_tcm_error];



    % Propagate to DV 2
    x_dv2_minus = stateProp(x_tcm_3_plus_r, dv2_time - t_tcm_3, simparams);

    % Compute corrected DV2
    [nom_dv_2_corrected, ~, exitflag, output] = fsolve(  @(nom_dv)fsolve_tcm_target_function2(nom_dv, x_dv2_minus, x(1:6,dv3_seg), dv3_time - dv2_time, simparams), nom_dv_2, options  );
    saveexitflag5(i) = exitflag;
    if exitflag ~=3
        badidx(i) = true;
        
        ppp=1;
    end
    dv2_correction_i(:,i) = nom_dv_2_corrected - nom_dv_2;

    dv2_corr_norm_diff_i(i) = norm(nom_dv_2_corrected) - norm(nom_dv_2);

    x_dv2_plus = x_dv2_minus + [zeros(3,1); nom_dv_2_corrected + randn(3,1) .* simparams.sig_dv_error];


    % Save the position at dv2 to compare with covariance ellipsoid
    dv2_position_i(:,i) = x_dv2_minus(1:3);




end



%% Follow-up calcs

dv2_pos = x(1:3,dv2_seg);
dv1_pos = x(1:3,2);


P_tcm_dv1 = T_dv2_dv1 * P_i_minus(:,:,dv1_event_idx) * T_dv2_dv1';
tcm_rss_dv1 = sqrt(trace(P_tcm_dv1));

P_tcm_dv2 = T_dv3_dv2 * P_i_minus(:,:,dv2_event_idx) * T_dv3_dv2';
tcm_rss_dv2 = sqrt(trace(P_tcm_dv2));


dv1_norm_i_mps = vecnorm(dv1_correction_i) * ndVel2kms * 1000;
dv1_std = std(dv1_norm_i_mps);

dv2_norm_i_mps = vecnorm(dv2_correction_i) * ndVel2kms * 1000;
dv2_std = std(dv2_norm_i_mps);


deltaV1 = deltaVs_nom0(:, 1 );

i_dv1 = deltaV1 / vecnorm(deltaV1);
dv1_correction_estimate = sqrt(i_dv1' * P_tcm_dv1 * i_dv1);
dv1_corr_std = std(dv1_corr_norm_diff_i);





% dv2_corr_norm_diff_i = rmoutliers(dv2_corr_norm_diff_i);


deltaV2 = deltaVs_nom0(:, 2 );

i_dv2 = deltaV2 / vecnorm(deltaV2);
dv2_correction_estimate = sqrt(i_dv2' * P_tcm_dv2 * i_dv2);


% dv2_dvi_std = std(dv2_corr_norm_diff_i) * ndVel2kms * 1000

dv2_dvi_std = std(dv2_corr_norm_diff_i(~[dv2_corr_norm_diff_i> 15*mean(dv2_corr_norm_diff_i)]))*ndVel2kms*1000;

dv2_dvi_mean = mean(dv2_corr_norm_diff_i(~[dv2_corr_norm_diff_i> 15*mean(dv2_corr_norm_diff_i)]))*ndVel2kms*1000;




%% figures



% change to dv1_corr_norm_diff_i

figure
histogram(dv1_corr_norm_diff_i*ndVel2kms*1000,'Normalization','probability')
hold on
% xline(tcm_rss_dv1*ndVel2kms*1000,'LineWidth',2,'Color','Black')
xline(dv1_correction_estimate*ndVel2kms*1000,'LineWidth',2,'Color','Red')
xline(dv1_corr_std*ndVel2kms*1000,'LineWidth',2,'Color','Blue')
xline(3*dv1_corr_std*ndVel2kms*1000,'LineWidth',2,'Color','Cyan')

legend('','i_d_v correction estimate','1\sigma of MC','3\sigma of MC')
xlabel('Correction \delta V (m/s)')
ylabel('Fraction in bin')
title('Corrected i \Delta V 1 Statistics')

% Histogram plot - corrected Delta V 1
figure
histogram(dv1_norm_i_mps,'Normalization','probability')
hold on
xline(tcm_rss_dv1*ndVel2kms*1000,'LineWidth',2,'Color','Black')
xline(dv1_correction_estimate*ndVel2kms*1000,'LineWidth',2,'Color','Red')
xline(dv1_std,'LineWidth',2,'Color','Blue')
xline(3*dv1_std,'LineWidth',2,'Color','Cyan')

legend('','P_t_c_m rss','i_d_v correction estimate','1\sigma of MC','3\sigma of MC')
xlabel('Correction \delta V (m/s)')
ylabel('Fraction in bin')
title('Corrected \Delta V 1 Statistics')




% Histogram plot - corrected Delta V 2
figure
histogram(dv2_norm_i_mps,'Normalization','probability')
hold on
xline(tcm_rss_dv2*ndVel2kms*1000,'LineWidth',2,'Color','Black')
xline(dv2_correction_estimate*ndVel2kms*1000,'LineWidth',2,'Color','Red')
xline(dv2_std,'LineWidth',2,'Color','Blue')
xline(3*dv2_std,'LineWidth',2,'Color','Cyan')

legend('','P_t_c_m rss','i_d_v correction estimate','1\sigma of MC','3\sigma of MC')
xlabel('Magnitude of difference in \delta V vectors (m/s)')
ylabel('Fraction in bin')
title('\Delta V 2 TCM Statistics')



% Histogram - dvi 2 corrected
figure
histogram(dv2_corr_norm_diff_i * ndVel2kms * 1000,'Normalization','probability')
hold on
% xline(tcm_rss_dv2*ndVel2kms*1000,'LineWidth',2,'Color','Black')
xline(dv2_correction_estimate*ndVel2kms*1000,'LineWidth',2,'Color','Red')
xline(dv2_dvi_std,'LineWidth',2,'Color','Blue')
xline(3*dv2_dvi_std,'LineWidth',2,'Color','Cyan')

xline(mean(dv2_corr_norm_diff_i)*ndVel2kms*1000,'LineWidth',2,'Color','Black')


% legend('','P_t_c_m rss','i_d_v correction estimate','1\sigma of MC','3\sigma of MC')
legend('','i_d_v correction estimate','1\sigma of MC','3\sigma of MC')
xlabel('Difference of magnitude of \Delta V vectors (m/s)')
ylabel('Fraction in bin')
title('\Delta V 2 Correction Statistics')






% histograms on components of corrections
figure
histogram(dv2_correction_i(1,:) * ndVel2kms * 1000,50,'Normalization','probability')


figure
histogram(dv2_correction_i(2,:) * ndVel2kms * 1000,50,'Normalization','probability')


figure
histogram(dv2_correction_i(3,:) * ndVel2kms * 1000,50,'Normalization','probability')



figure;
if simparams.perform_correction    
%     [~,tcm_idx0] = opt_multiple_tcm(simparams.x0, deltaVs_nom0, t0, t_s0, stm_t0, stm_t_i0, simparams);
    plotMultiSegTraj(simparams.x0, x_t0, t_s0, simparams, tcm_idx0);
else
    plotMultiSegTraj(simparams.x0, x_t0, t_s0, simparams);
end
axis equal;
% Adding position at DV2 to trajectory figure
hold on;
xlim([1-2*mu, 1])
zlim([-mu mu])
ylim([-mu mu])
plot3(dv2_position_i(1,:),dv2_position_i(2,:),dv2_position_i(3,:),'.','MarkerSize',25)
view([115 10])



f1 = gca;
f2 = figure;
a2 = copyobj(f1,f2);
hold on;
xlim([dv2_pos(1) - .00001, dv2_pos(1) + .00001])
zlim([dv2_pos(3) - .00001, dv2_pos(3) + .00001])
ylim([dv2_pos(2) - .00001, dv2_pos(2) + .00001])

% for i = 1:nRuns
%     drawVector(dv2_pos', dv2_correction_i(:,i)', 'green');
%     drawVector(dv2_position_i(:,i)', dv2_correction_i(:,i)' * .01, 'green',1);
    
% end


% pause;

% Also adding the dispersion covariance ellipsoid to the figure
[ex, ey, ez] = ellipsoid(P_i_minus(1:3,1:3,dv2_event_idx), 100);
surf(dv2_pos(1)+ex*3, dv2_pos(2)+ey*3, dv2_pos(3)+ez*3, 'FaceAlpha',.5,'EdgeAlpha',0);

% pause;

% Adding the TCM covariance ellipsoid to the figure
[ex, ey, ez] = ellipsoid(P_tcm_dv2, 100);
surf(dv2_pos(1)+ex*.005, dv2_pos(2)+ey*.005, dv2_pos(3)+ez*.005, 'FaceAlpha',.5,'EdgeAlpha',0);


% pause;
% Add the nominal DV unit vector
drawVector(dv2_pos', i_dv2'*.00015,'blue');
% pause;

drawVector(dv2_position_i', dv2_correction_i' * .005, 'green',1);









% Histogram plot - TCM 1
tcm1_norm_i = vecnorm(tcm_1_corrected_i) * ndVel2kms * 1000;
P_tcm1 = T1 * P_i_minus(:,:,2) * T1';
tcm1_rss = sqrt(trace(P_tcm1));
tcm1_997 = prctile(tcm1_norm_i,99.7);

figure;
histogram(tcm1_norm_i,'Normalization','probability')

xline(3*tcm1_rss*ndVel2kms*1000,'LineWidth',2,'Color','Black')
xline(tcm1_997,'LineWidth',2,'Color','Blue')

legend('','3P_t_c_m rss','99.7% of MC')
xlabel('TCM 1 \delta V (m/s)')
ylabel('Fraction in bin')
title('TCM 1 Statistics')




% Histogram plot - TCM 2
tcm2_norm_i = vecnorm(tcm_2_corrected_i) * ndVel2kms * 1000;
P_tcm2 = T2 * P_i_minus(:,:,3) * T2';
tcm2_rss = sqrt(trace(P_tcm2));
tcm2_997 = prctile(tcm2_norm_i,99.7);

figure;
histogram(tcm2_norm_i,'Normalization','probability')

xline(3*tcm2_rss*ndVel2kms*1000,'LineWidth',2,'Color','Black')
xline(tcm2_997,'LineWidth',2,'Color','Blue')

legend('','3P_t_c_m rss','99.7% of MC')
xlabel('TCM 2 \delta V (m/s)')
ylabel('Fraction in bin')
title('TCM 2 Statistics')





% Histogram plot - TCM 3
tcm3_norm_i = vecnorm(tcm_3_corrected_i) * ndVel2kms * 1000;
P_tcm3 = T3 * P_i_minus(:,:,4) * T3';
tcm3_rss = sqrt(trace(P_tcm3));
tcm3_997 = prctile(tcm3_norm_i,99.7);

figure;
histogram(tcm3_norm_i,'Normalization','probability')

xline(3*tcm3_rss*ndVel2kms*1000,'LineWidth',2,'Color','Black')
xline(tcm3_997,'LineWidth',2,'Color','Blue')

legend('','3P_t_c_m rss','99.7% of MC')
xlabel('TCM 3 \delta V (m/s)')
ylabel('Fraction in bin')
title('TCM 3 Statistics')






%% Fmincon call via output function
% tic
% [x_opt,J,history,searchdir,exitflag,output] = runfmincon(simparams);
% toc


%% View optimal solution

% % Calculate optimal TCM time and delta V values
% [stm_i, stt_i, x_i_f, x_t, stm_t, stt_t_i, t, t_s, stm_t_i]  = createStateStmSttHistory(x_opt, simparams);
% % Calculate total impulsive delta V 
% [deltaV, deltaVs_nom] = calcDeltaV(x_opt, x_i_f, simparams);
% [tcm_time,tcm_idx,min_tcm_dv,~,~,tcm_dv_each] = opt_multiple_tcm(x_opt, deltaVs_nom, t, t_s, stm_t, stm_t_i, simparams);
% totalDV = deltaV + 3*min_tcm_dv
% 
% 
% 
% solfig = figure;
% plotMultiSegTraj(x_opt, x_t, t_s, simparams, tcm_idx);
% % plotMultiSegTraj(x_opt, x_t, t_s, simparams);
% title('optim_solution')
% solfig.CurrentAxes.Title.Visible="off";





%% Save outputs

% if output.firstorderopt < 4 && exitflag ~= -2 && saveOutput
% if exitflag ~= -2 && saveOutput
if saveOutput


    formatOut = 'yyyymmdd_HHMM.SS';
    dateString = datestr(now,formatOut);
    outputPath = strcat('../sims/',dateString,'_',savename);
    mkdir(outputPath);
    outputPathName = strcat('./',outputPath,'/video.avi');
 
    % Save workspace    
    saveallfigs(outputPath,0,1);
    save(strcat('./',outputPath,'/workspace.mat'));
    
    
    % save main file backup
    FileNameAndLocation=[mfilename('fullpath')];
    newbackup=strcat(outputPath,'/main_mc_backup.m');
    % newbackup=sprintf('%sbackup.m',FileNameAndLocation);
    currentfile=strcat(FileNameAndLocation, '.m');
    copyfile(currentfile,newbackup);
    
       
    % save simulation parameters backup
    newbackup=strcat(outputPath,'/initialize_simparams_backup.m');
    copyfile(strcat(init_fn,'.m'),newbackup);


end








% Save iteration history .fig only
% iterhist = figure;
% plotIterationHistory(history.x,simparams)
% title('iteration_history')
% iterhist.CurrentAxes.Title.Visible="off";
% 
% savefig(iterhist, strcat(outputPath,'/iteration_history.fig'));


%% debug
% [stm_i0, stt_i0, x_i_f0, x_t0, stm_t0, stt_t_i, t0, t_s0, stm_t_i0]  = createStateStmSttHistory(x, simparams);
% % [tcm_time0, tcm_idx0, min_tcm_dv0, ~, ~, tcm_dv_each0] = opt_multiple_tcm(x, t0, t_s0, stm_t0, simparams); % inputs: x, t, t_s, stm_t, stm_t_i, simparams
% figure
% plotMultiSegTraj(x, x_t0, t_s0, simparams);







%% investigating apse constraint
% v_I = zeros(3,length(t));
% 
% rdotv = zeros(1,length(t));
% r_m_sc_mag = rdotv;
% 
% for i = 1:length(t)
% 
%     omega_SI_nd = [0; 0; 1];
%     r_m = [1-mu; 0; 0];
% 
%     % Get inertial velocity
%     
%     v_S = x_t(i,4:6)';
%     r_S = x_t(i,1:3)';
%     v_I(:,i) = v_S + cross(omega_SI_nd, r_S);
% 
% 
%     r_m_sc = r_S - r_m;
%     r_m_sc_mag(i) = vecnorm(r_m_sc);
%     i_r_m_sc = r_m_sc/vecnorm(r_m_sc);
% 
%     i_v_I = v_I(:,i)/vecnorm(v_I(:,i));
% 
%     rdotv(i) = i_r_m_sc' * i_v_I;
% 
%     
% end
% 
% 
% figure
% plot(rdotv)
% hold on;
% yline(0)
% plot(r_m_sc_mag)


% 
% test_dot = abs(rdotv(4100:4300));
% [minv,minidx] = min(test_dot);
% 
% rdotv(4100+minidx-1)
% 
% 
% 
% plot3(x_t(4100+minidx-1,1),x_t(4100+minidx-1,2),x_t(4100+minidx-1,3),'.','MarkerSize',25)


