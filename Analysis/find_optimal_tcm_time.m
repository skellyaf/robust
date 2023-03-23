clc; close all; clear;

%% load a LEO coast to Hohmann xfer to GEO 
load('C:\Users\skell\OneDrive - USU\Documents\code_repos\multiseg_opt\sims\20220706_1627.55_leo2leo_smallP_objandconstraintGrad_testsnoSTT_EXAMPLE1\workspace.mat')
addpath(genpath('../'))
% Traj is in x_opt
x = x_opt;

%% Define P_initial
% Initial uncertainty
% Small
% sig_pos = 10 / 1e3; % Position +/- 10 m in all 3 direction
% sig_vel = .01 / 1e3; % Velocity +/- 1 cm/s in all 3 directions
% Medium
% sig_pos = 1000 / 1e3; % Position +/- 1 km in all 3 direction
% sig_vel = 1 / 1e3; % Velocity +/- 1 m/s in all 3 directions
% Large
% sig_pos = 10000 / 1e3; % Position +/- 10 km in all 3 direction
% sig_vel = 10 / 1e3; % Velocity +/- 10 m/s in all 3 directions
% Huge
% sig_pos = 100000 / 1e3; % Position +/- 100 km in all 3 direction
% sig_vel = 100 / 1e3; % Velocity +/- 100 m/s in all 3 directions


sig_pos = 1; % Position +/- 10 m in all 3 direction
sig_vel = 0; % Velocity +/- 1 cm/s in all 3 directions
% P_initial = diag([sig_pos^2 sig_pos^2 sig_pos^2 sig_vel^2 sig_vel^2 sig_vel^2]);
P_initial = diag([0 sig_pos^2 sig_pos^2  sig_vel^2 sig_vel^2 sig_vel^2]);

%% Parse segments
% # of parameters per segment
m = 7;

% Calculate the number of corrections
% This works currently as long as there aren't m or more corrections
num_corrections = mod(length(x),m);
n = ( length(x) - num_corrections ) / m;

% Array of corrective maneuver times that is num_corrections long
t_corr = x(end-num_corrections+1:end);

% Reshaping x so each column is a segment initial state and duration
x = reshape(x(1:end-num_corrections),m,n);


%% Propagate entire trajectory, save dynamics at each time step

for i = 1:n

    x_i_initial = x(1:6,i);
    delta_t = x(7,i);

    if delta_t > 0
        [x_i_final, stm_i(:,:,i), xstm_t_i, t_i] = stateStmProp(x_i_initial, delta_t, mu);
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
    dstmNCdtc = -stmNC * r2bp_A_matrix(x_c);
    dTdtc = T_partial(stmNC,dstmNCdtc);
    dstmC0dtc = r2bp_A_matrix(x_c) * stmC0;
    corrGrad(i) = sigDvPartial(T,stmC0,P_initial,dTdtc,dstmC0dtc);   
    

end






zci = @(v) find(diff(sign(v)));

ZeroX = @(x0,y0,x1,y1) x0 - (y0.*(x0 - x1))./(y0 - y1);

zxidx = zci(corrGrad);

t_zerox = [];
corrMag = [];


for m = 1:length(zxidx)
    idx1 = zxidx(m);
    idx2 = idx1+1;
    t_1 = t(idx1);
    t_2 = t(idx2);

    y1 = corrGrad(idx1);
    y2 = corrGrad(idx2);

    % Are the slopes the same sign or different signs before, across, and after the
    % correction? (ie is the sign change a zero crossing or asymptote?)

    slope1 = corrGrad(idx1) - corrGrad(idx1-1);
    slope2 = corrGrad(idx2) - corrGrad(idx1);
    slope3 = corrGrad(idx2+1) - corrGrad(idx2);

    % If it is zero crossing
    if abs(  sign(slope1) + sign(slope2) + sign(slope3)   ) == 3
        t_zerox(end+1) = ZeroX(t_1,y1,t_2,y2);

        % Calculate correction magnitude
        x_currentCorrection = [x(:); t_zerox(end)];
        [corrMag(end+1)] = xfer_obj_corr_grad(x_currentCorrection,P_initial, mu);



    end


    
end

[minDV, minIdx] = min(corrMag);
t_c_opt = t_zerox(minIdx)


plot(t, corrGrad)
hold on
yline(0)
limit = trimmean(corrGrad,20);
ylim([-limit limit])

plot(t_zerox, 0*t_zerox,'.','MarkerSize',25)

plot(t_c_opt,0,'o','MarkerSize',10)


% 
% corrMag = [];
% for m = 1:length(t_zerox)
%     x(end) = t_zerox(m);
%     [corrMag(end+1)] = xfer_obj_corr_grad(x,P_initial, mu)
%     if abs(corrGrad(end)) > 100
%         corrMag(end) = NaN;
%     end
% end
% 
% 
% Ztcm(i,j) = min(corrMag);



