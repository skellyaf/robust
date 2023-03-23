function [tcm_gradient] = calc_tcm_gradient(x, x_i_f, stm_i, stt_i, stm_t, stt_t_i, t, t_s, tcm_time, dvR3sigma_tr, dvV3sigma_tr, simparams)
%calc_tcm_gradient Computes and returns the analytical tcm gradient

%   Detailed explanation goes here


m = simparams.m;
n = simparams.n;
mu = simparams.mu;
P_initial = simparams.P_initial;

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




maneuverSegments = simparams.maneuverSegments;
num_corrections = 1;

% One gradient vector for each corrective maneuver, to be summed later
corrgradJ = zeros(n*m+num_corrections,num_corrections);

% Gradient of the TCM portion of the objective function
tcm_gradient = zeros(n*m, 1);

% Test which segment has the correction
corrSeg = t_s(t==tcm_time);


% Loop over each segment, calculate gradients of tcm wrt x and delta T
% for i = 1:n
% 
%     % dTCM / dX0,i
% 
% 
%     % dTCM / d Delta t,i
% 
% end

stmC0 = stm_t(:,:,t==tcm_time);
if corrSeg > 1
    t_corrSegStart = sum(x(7,1:corrSeg-1));
else
    t_corrSegStart = 0;
end

t_corrSegEnd = sum(x(7,1:corrSeg));

stmCs0 = stm_t(:,:,t==t_corrSegStart);
stm0Cs = -J*stmCs0'*J;
stmCCs = stmC0 * stm0Cs;

stm_leg_1 = stmCCs;

% stmCsC = -J * stmCCs' * J;

stmCf0 = stm_t(:,:,t==t_corrSegEnd);
stm0C = -J * stmC0' * J;
stmCfC = stmCf0 * stm0C;
stm_leg_2 = stmCfC;

stmNC = stmN0 * stm0C;
T = [-inv( stmNC(1:3,4:6) ) * stmNC(1:3,1:3), -eye(3)];



% Calculate the STT from the beginning of the correction segment to the
% correction
stt_t_corrSeg = stt_t_i{corrSeg}(:,:,:,2:end);
t_corrSeg = t(t_s == corrSeg);

stt_leg_1 = stt_t_corrSeg(:,:,:,t_corrSeg == tcm_time);





for j = 1:length(corrSeg) % Placeholder for more than 1 correction

    for i = 1:n

        if i > 1
            stmIm0 = stmCombine(stm_i, 1, i-1);
        else
            stmIm0 = eye(6);
        end

        if i < corrSeg(j)
            %corrgradJ(:,j)



            if i + 1 < corrSeg(j)
                stmCmIp = stmCombine(stm_i, i + 1, corrSeg(j) - 1);
                sttCmIp = sttCombine(stm_i, stt_i, i + 1, corrSeg(j) - 1);
                sttCmIp = tensorCombine(stmCmIp, sttCmIp, stm_leg_1, stt_leg_1);
            else
                stmCmIp = eye(6);
                sttCmIp = stt_leg_1;
            end

            % Partial of variance wrt state x
            dstmC0dx = tmult(stm_leg_1 * stmCmIp, tmult(stt_i(:,:,:,i), stmIm0, [0 0]), [0 0]);
            corrgradJ((i-1)*7+1:(i-1)*7+6,j) = sigDvPartial(T,stmC0,P_initial,zeros(3,6,6),dstmC0dx);

            % Partial of variance wrt segment time, delta t
%                     dstmC0ddt = stm_leg_1 * stmCmIp * r2bp_A_matrix(x_i_f(:,i)) * stm_i(:,:,i) * stmIm0 + tensorVectorMult(sttCmIp, r2bp_de(1, x(1:6,i+1))) * stm_i(:,:,i) * stmIm0;
            dstmC0ddt = stm_leg_1 * stmCmIp * r2bp_A_matrix(x_i_f(:,i), mu) * stm_i(:,:,i) * stmIm0;
%                     dstmNCdti = stmNCp * tensorVectorMult(stt_leg_2, r2bp_de(1, x_c)) * stm_leg_1 * stm_i(:,:,i)
%                     dstmNCdti = stmNC * stm_leg_1 * stmCmIp * tensorVectorMult(stt_i(:,:,:,i), r2bp_de(1, x(1:6,i)))
%                     dstmNCdti = stmNC * stm_leg_1 * stmCmIp * r2bp_A_matrix(x_i_f(:,i))
%                     dTddt = T_partial(stmNC, dstmNCdti)


%                     corrgradJ((i-1)*7+7,j) = sigDvPartial(T,stmC0,P_initial,dTddt,dstmC0ddt)
            corrgradJ((i-1)*7+7,j) = sigDvPartial(T,stmC0,P_initial,zeros(3,6),dstmC0ddt);
        elseif i == corrSeg(j)
            % When i = the correction segment


            dstmNCdxi_num = zeros(6,6,6);
            dx = 1e-8;
            for k = 1:6
                dxi = zeros(6,1);    
                dxi(k) = dx;
                dx_c_initial = x_c_initial + dxi; % varying x_initial

                [dx_c, dstm_leg_1] = stateStmProp(dx_c_initial, leg_t_c, simparams);
                [~, dstm_leg_2] = stateStmProp(dx_c, delta_t_c - leg_t_c, simparams);
                
                dstmNCdxi_num(:,:,k) = (stmNCp*dstm_leg_2 - stmNCp*stm_leg_2)./dx;
            end


            % Can't figure this one out - doing it numerically
            % above for now
            dstmNCdxi = tmult(stmNC, stt_leg_1);





            dTdxi = T_partial(stmNC, dstmNCdxi_num);
            dstmC0dxi = tmult(stt_leg_1,stmIm0);
            corrgradJ((i-1)*7+1:(i-1)*7+6,j) = sigDvPartial(T,stmC0,P_initial,dTdxi,dstmC0dxi);


            % WRT time
            if i == n
                ppp=1;
            end
%                     if i < n
            if 1
%                         dstmNCddti2 = stmNCp * r2bp_A_matrix(x(:,i+1)) * stm_leg_2;
                % Investigate the difference between above and
                % below
                dstmNCddti = stmNCp * r2bp_A_matrix(x_i_f(:,i), mu) * stm_leg_2;
            else
                assert(0,'Need to do work here');
            end

            dTddti = T_partial(stmNC, dstmNCddti);


            corrgradJ((i-1)*7+7,j) = sigDvPartial(T,stmC0,P_initial,dTddti,zeros(6,6));
        else
            % When i > correction segment
            if i + 1 <= n
                stmNIp = stmCombine(stm_i,i+1,n);
%                         stmIC = stmCombine(stm_i,corrSeg(j),i+1)*stm_leg_2;
            else
                stmNIp = eye(6);
%                         stmIC = stm_leg_2;
            end



            if i-1 > corrSeg(j)
                stmImC = stmCombine(stm_i,corrSeg(j)+1,i-1) * stm_leg_2;
            else
                stmImC = stm_leg_2;
            end

            stmIC = stm_i(:,:,i)*stmImC;
            

            dstmNCdxi = tmult(stmNIp, tmult(stt_i(:,:,:,i), stmImC));
            dTdxi = T_partial(stmNC, dstmNCdxi);
            corrgradJ((i-1)*7+1:(i-1)*7+6,j) = sigDvPartial(T,stmC0,P_initial,dTdxi,zeros(6,6,6));

              
            
            dstmNCddti = stmNIp * r2bp_A_matrix(x_i_f(:,i), mu) * stmIC;
            dTddti = T_partial(stmNC, dstmNCddti);
            corrgradJ((i-1)*7+7,j) = sigDvPartial(T,stmC0,P_initial,dTddti,zeros(6,6));

        end
    end
end

% GRADIENT OF CORRECTION TIME ONLY
if testvar == 0
    ppp=1;
end

dstmNCdtc = -stmNC * r2bp_A_matrix(x_c, mu);
dTdtc = T_partial(stmNC,dstmNCdtc);
dstmC0dtc = r2bp_A_matrix(x_c, mu) * stmC0;
corrgradJ(end,1) = sigDvPartial(T,stmC0,P_initial,dTdtc,dstmC0dtc);   





tcm_gradient = 9*1/2*dvR3sigma_tr^(-1)*corrgradJ;

end