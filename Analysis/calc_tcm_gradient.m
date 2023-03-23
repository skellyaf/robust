function [tcm_gradient] = calc_tcm_gradient(x, x_t, x_i_f, stm_i, stt_i, stm_t, stt_t_i, t, t_s, tcm_time, dvR3sigma_tr, dvV3sigma_tr, simparams)
%calc_tcm_gradient Computes and returns the analytical tcm gradient

%   Detailed explanation goes here

% Variable extraction
m = simparams.m;
n = simparams.n;
mu = simparams.mu;
P_initial = simparams.P_initial;
maneuverSegments = simparams.maneuverSegments;
num_corrections = 1;

% Reshaping x so each column is a segment initial state and duration
x = reshape(x,m,n);

% Symplectic unit matrix
J = [zeros(3,3), eye(3,3); -eye(3,3), zeros(3,3)];

% STM from beginning of trajectory to the state being targeted
if simparams.target_final_maneuver
    final_seg = simparams.maneuverSegments(end);
    target_time = sum(x(7,1:final_seg - 1));
    target_idx = find(target_time == t) ;
    stmN0 = stm_t(:,:,target_idx);
else
    stmN0 = stm_t(:,:,end);
end

% One gradient vector for each corrective maneuver, to be summed later
corrgradJ = zeros(n*m,num_corrections); % modified, t_c is no longer a parameter
%%%%% replacing corrgradJ with tcm_gradient_r and tcm_gradient_v

% Gradient of the TCM portion of the objective function
tcm_gradient_r = zeros(n*m, num_corrections);
tcm_gradient_v = zeros(n*m, num_corrections);

%% Setup for logic choices and STM portions
% Test which segment has the correction
corrSeg = t_s(t==tcm_time);

% Get the nominal state at the correction
x_tcm = x_t(t==tcm_time,:)';

% Extract STMs for calculations
% STM from the beginning of entire traj (0) to the correction (C) = stmC0
stmC0 = stm_t(:,:,t==tcm_time);

% Get the time the correction segment starts = t_corrSegStart
if corrSeg > 1
    t_corrSegStart = sum(x(7,1:corrSeg-1));
else
    t_corrSegStart = 0;
end

% Get the time the correction segment ends = t_corrSegEnd
t_corrSegEnd = sum(x(7,1:corrSeg));

% Get the STM from the beginning of the trajectory to the beginning of the 
% correction segment = stmCs0
stmCs0 = stm_t(:,:,t==t_corrSegStart);
stm0Cs = -J*stmCs0'*J;

% STM from the start of the correction segment to the correction = stmCCs
stmCCs = stmC0 * stm0Cs;
stm_leg_1 = stmCCs;
stmCsC = -J * stmCCs' * J;

% Get the STM from the correction to the end of the correction seg = stmCfC
stmCf0 = stm_t(:,:,t==t_corrSegEnd);
stm0C = -J * stmC0' * J;
stmCfC = stmCf0 * stm0C;
stm_leg_2 = stmCfC;

% Calculate the STM from the correction to the end of the traj = stmNC
% and the T matrix (used for TCMr calculation)
stmNC = stmN0 * stm0C;
T = [-inv( stmNC(1:3,4:6) ) * stmNC(1:3,1:3), -eye(3)];

% Calculate the W and L matrices (used for TCMv calculation)
W = [stmNC(4:6,4:6), -stmNC(4:6,1:3)];
L = [W*T', zeros(3,3)];

% Calculate the STT from the beginning of the correction segment to the
% correction
stt_t_corrSeg = stt_t_i{corrSeg}(:,:,:,2:end);
t_corrSeg = t(t_s == corrSeg);

stt_leg_1 = stt_t_corrSeg(:,:,:,t_corrSeg == tcm_time);

% Calculate the STM from the beginning of the segment after the correction 
% to the target
% stmN0 already exists
% stmCf0 already exists
stm0Cf = -J * stmCf0' * J;
stmNCf = stmN0 * stm0Cf;

% Pre-calculate the four leading terms inside the variance partial eqn:
%%% dTCMr/dx0i = 2* (  tr(T*stmc0*P0*stmC0 * dT/dx) + tr(T'*T*stmC0*P0 * sttC0)  )
% In the 





%% Loop through the segments / calculate partial derivatives / gradients

for j = 1:length(corrSeg) % Placeholder for more than 1 correction

    for i = 1:n

        % FIND:
        
        % -- dTdxi
        % -- dWdxi
        % -- dLdxi
        % ---- The above 3 need the following, which vary by segment
        % ---- dstmC0dxi
        % ---- dstmFCdxi

        % Get the STT for the i-th segment
        stti = stt_i(:,:,:,i);

        t_i = t(t_s==i); % Time elements corresponding to seg i

        % Logic case 1: segment i is before the segment containing the TCM
        if i < corrSeg(j)
            % ---- dstmC0dxi
            % Get the STM before segment i (if it exists)
            stmi0 = stm_t(:,:,t==t_corrSegStart);

            % Get the STM from the end of the i-th segment to the beginning
            % of the correction segment: stmCsif
            
            stmif0 = stm_t(:,:,t==t_i(end));
            stm0if = -J * stmif0' * J;
            stmCsif = stmCs0 * stm0if;
            stmCif = stmCCs * stmCsif;

            dstmC0dxi = tmult(stmCif, tmult(stti, stmi0));

            % ---- dstmC0ddti
            dstmC0ddti = stmCif * r2bp_A_matrix(x_i_f(:,i), mu) * stmif0;

            % ---- dstmNCdxi
            dstmNCdxi = zeros(6,6,6);

            % ---- dstmNCddti
            dstmNCddti = zeros(6,6);


        % Logic case 2: segment i contains the TCM
        elseif i == corrSeg(j)
            % ---- dstmNCdxi
            % Need to calculate the ISTT of sttCCs
            % First calculate sttCCs (from cell structure stt_t_i)
            stt_t_C = stt_t_i{i}(:,:,:,2:end);
            sttCCs = stt_t_C(:,:,:,t_i==tcm_time);
            % Invert sttCCs to get sttCsC
            sttCsC = invert_stt(stmCCs, sttCCs);

            dstmNCdxi = tmult(stmNCf, tmult(stti,stmCsC) + tmult(stm_i(:,:,i), sttCsC));

            % ---- dstmNCddti
            % Need state at end of correction segment
            x_c_f = x_i_f(:,i);
            dstmNCddti = stmNCf * r2bp_A_matrix(x_c_f, mu) * stmCfC;

            % ---- dstmC0dxi
            dstmC0dxi = tmult(sttCCs,stmCs0);

            % ---- dstmC0ddti
            dstmC0ddti = r2bp_A_matrix(x_tcm, mu) * stmC0;
            dstmC0ddti = zeros(6,6);


        % Logic case 3: segment i is after the segment containing the TCM
        % but not after the target
        elseif target_time > sum(x(7,1:i-1))

            % ---- dstmNCdxi
            % Get STM from the end of i-th segment to the target
            % Already have stmN0, get stmif0

            stmif0 = stm_t(:,:,t==t_i(end));
            stm0if = -J * stmif0' * J;
            stmNif = stmN0 * stm0if;

            % Already have STT from beginning to end of segment i, stti

            % Get STM from the correction to the beginning of segment i
            % (end of segment i-1)
            % Already have stmC0 and stm0C
            % Get stmi0
            t_i_minus_1 = t(t_s==i-1);
            stmi0 = stm_t(:,:,t==t_i_minus_1(end));
            stmiC = stmi0 * stm0C;

            dstmNCdxi = tmult(stmNif, tmult(stti, stmiC));

            % ---- dstmNCddti
            dstmNCddti = stmNif * r2bp_A_matrix(x_i_f(:,i), mu) * stm_i(:,:,i) * stmiC;

            % ---- dstmC0dxi
            dstmC0dxi = zeros(6,6,6);

            % ---- dstmC0ddti
            dstmC0ddti = zeros(6,6);

        % Logic case 4: segment i is after the target (coast segment)
        else 
            % ---- dstmNCdxi
            dstmNCdxi = zeros(6,6,6);
            % ---- dstmNCddti
            dstmNCddti = zeros(6,6);
            % ---- dstmC0dxi
            dstmC0dxi = zeros(6,6,6);
            % ---- dstmC0ddti
            dstmC0ddti = zeros(6,6);

        end


        % Assemble TCMr gradient with calculated pieces above
        dTdxi = T_partial(stmNC, dstmNCdxi);
        tcm_gradient_r((i-1)*7+1:(i-1)*7+6,j) = sigDvPartial(T,stmC0,P_initial,dTdxi,dstmC0dxi);

        dTddti = T_partial(stmNC, dstmNCddti);
        tcm_gradient_r((i-1)*7+7,j) = sigDvPartial(T,stmC0,P_initial,dTddti,dstmC0ddti);

        % Assemble TCMv gradient with calculated pieces above
        dLdxi = L_partial(W, T, dstmNCdxi, dTdxi);
        tcm_gradient_v((i-1)*7+1:(i-1)*7+6,j) = sigDvPartial(L,stmC0,P_initial,dLdxi,dstmC0dxi);

        dLddti = L_partial(W, T, dstmNCddti, dTddti);
        tcm_gradient_v((i-1)*7+7,j) = sigDvPartial(L,stmC0,P_initial,dLddti,dstmC0ddti);





    end

end

% assert(0)
% for j=1
%     for i=1:n
%         x_c_initial = x(1:6,corrSeg);
%         if i > 1
%             stmIm0 = stmCombine(stm_i, 1, i-1);
%         else
%             stmIm0 = eye(6);
%         end
% 
% 
%         % Logic case 1: segment i is before the segment containing the TCM
%         if i < corrSeg(j)
% 
%             
% 
% 
%             % Already have STM from beginning of traj to correction: stmC0
%             % Sensitivity 
%             
%             if i + 1 < corrSeg(j)
%                 stmCmIp = stmCombine(stm_i, i + 1, corrSeg(j) - 1);
%                 sttCmIp = sttCombine(stm_i, stt_i, i + 1, corrSeg(j) - 1);
%                 sttCmIp = tensorCombine(stmCmIp, sttCmIp, stm_leg_1, stt_leg_1);
%             else
%                 stmCmIp = eye(6);
%                 sttCmIp = stt_leg_1;
%             end
% 
%             % Partial of variance wrt state x
%             dstmC0dx = tmult(stm_leg_1 * stmCmIp, tmult(stt_i(:,:,:,i), stmIm0, [0 0]), [0 0]);
%             dT = zeros(3,6,6);
%             corrgradJ((i-1)*7+1:(i-1)*7+6,j) = sigDvPartial(T,stmC0,P_initial,dT,dstmC0dx);
% 
% 
% 
%             % Partial of variance wrt segment time, delta t
% %                     dstmC0ddt = stm_leg_1 * stmCmIp * r2bp_A_matrix(x_i_f(:,i)) * stm_i(:,:,i) * stmIm0 + tensorVectorMult(sttCmIp, r2bp_de(1, x(1:6,i+1))) * stm_i(:,:,i) * stmIm0;
%             dstmC0ddt = stmC0 * r2bp_A_matrix(x_i_f(:,i), mu) * stm_i(:,:,i) * stmCs0;
% %                     dstmNCdti = stmNCp * tensorVectorMult(stt_leg_2, r2bp_de(1, x_c)) * stm_leg_1 * stm_i(:,:,i)
% %                     dstmNCdti = stmNC * stm_leg_1 * stmCmIp * tensorVectorMult(stt_i(:,:,:,i), r2bp_de(1, x(1:6,i)))
% %                     dstmNCdti = stmNC * stm_leg_1 * stmCmIp * r2bp_A_matrix(x_i_f(:,i))
% %                     dTddt = T_partial(stmNC, dstmNCdti)
% 
% 
% %                     corrgradJ((i-1)*7+7,j) = sigDvPartial(T,stmC0,P_initial,dTddt,dstmC0ddt)
%             corrgradJ((i-1)*7+7,j) = sigDvPartial(T,stmC0,P_initial,zeros(3,6),dstmC0ddt);
%         elseif i == corrSeg(j)
%             % When i = the correction segment
% 
% 
%             dstmNCdxi_num = zeros(6,6,6);
%             dx = 1e-8;
%             for k = 1:6
%                 dxi = zeros(6,1);    
%                 dxi(k) = dx;
%                 dx_c_initial = x_c_initial + dxi; % varying x_initial
% 
%                 [dx_c, ~] = stateStmProp(dx_c_initial, leg_t_c, simparams);
%                 [~, dstm_leg_2] = stateStmProp(dx_c, delta_t_c - leg_t_c, simparams);
%                 
%                 dstmNCdxi_num(:,:,k) = (stmNCp*dstm_leg_2 - stmNCp*stm_leg_2)./dx;
%             end
% 
% 
%             % Can't figure this one out - doing it numerically
%             % above for now
%             dstmNCdxi = tmult(stmNC, stt_leg_1);
% 
% 
% 
% 
% 
%             dTdxi = T_partial(stmNC, dstmNCdxi_num);
%             dstmC0dxi = tmult(stt_leg_1,stmCs0);
%             corrgradJ((i-1)*7+1:(i-1)*7+6,j) = sigDvPartial(T,stmC0,P_initial,dTdxi,dstmC0dxi);
% 
% 
%             % WRT time
%             if i == n
%                 ppp=1;
%             end
% %                     if i < n
%             if 1
% %                         dstmNCddti2 = stmNCp * r2bp_A_matrix(x(:,i+1)) * stm_leg_2;
%                 % Investigate the difference between above and
%                 % below
%                 dstmNCddti = stmNCp * r2bp_A_matrix(x_i_f(:,i), mu) * stm_leg_2;
%             else
%                 assert(0,'Need to do work here');
%             end
% 
%             dTddti = T_partial(stmNC, dstmNCddti);
% 
% 
%             corrgradJ((i-1)*7+7,j) = sigDvPartial(T,stmC0,P_initial,dTddti,zeros(6,6));
%         else
%             % When i > correction segment
%             if i + 1 <= n
%                 stmNIp = stmCombine(stm_i,i+1,n);
% %                         stmIC = stmCombine(stm_i,corrSeg(j),i+1)*stm_leg_2;
%             else
%                 stmNIp = eye(6);
% %                         stmIC = stm_leg_2;
%             end
% 
% 
% 
%             if i-1 > corrSeg(j)
%                 stmImC = stmCombine(stm_i,corrSeg(j)+1,i-1) * stm_leg_2;
%             else
%                 stmImC = stm_leg_2;
%             end
% 
%             stmIC = stm_i(:,:,i)*stmImC;
%             
% 
%             dstmNCdxi = tmult(stmNIp, tmult(stt_i(:,:,:,i), stmImC));
%             dTdxi = T_partial(stmNC, dstmNCdxi);
%             corrgradJ((i-1)*7+1:(i-1)*7+6,j) = sigDvPartial(T,stmC0,P_initial,dTdxi,zeros(6,6,6));
% 
%               
%             
%             dstmNCddti = stmNIp * r2bp_A_matrix(x_i_f(:,i), mu) * stmIC;
%             dTddti = T_partial(stmNC, dstmNCddti);
%             corrgradJ((i-1)*7+7,j) = sigDvPartial(T,stmC0,P_initial,dTddti,zeros(6,6));
% 
%         end
%     end
% end
% 
% % GRADIENT OF CORRECTION TIME ONLY
% if testvar == 0
%     ppp=1;
% end
% 
% dstmNCdtc = -stmNC * r2bp_A_matrix(x_c, mu);
% dTdtc = T_partial(stmNC,dstmNCdtc);
% dstmC0dtc = r2bp_A_matrix(x_c, mu) * stmC0;
% corrgradJ(end,1) = sigDvPartial(T,stmC0,P_initial,dTdtc,dstmC0dtc);   
% 




% tcm_gradient = 9*1/2*dvR3sigma_tr^(-1)*corrgradJ;







tcm_gradient = 9/2*(dvR3sigma_tr^(-1)*tcm_gradient_r + dvV3sigma_tr^(-1)*tcm_gradient_v);

end