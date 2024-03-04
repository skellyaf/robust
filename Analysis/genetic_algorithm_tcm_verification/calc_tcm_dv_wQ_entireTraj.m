function tcm_dv_total = calc_tcm_dv_wQ_entireTraj(tcm_idx, x, traj, deltaVs_nom, simparams)



tcm_time = traj.t(tcm_idx)';

[Q_k_km1] = calc_Q_events(traj, x, tcm_time, simparams);

[~, tcm_dv_total] = calc_covariance_wQ_tcmdv_v3(x(:), traj, tcm_time, 1, deltaVs_nom, simparams.P_initial, Q_k_km1, simparams);











%     t_eval = t(start_idx:target_idx);
%     t_s_eval = t_s(start_idx:target_idx);
%     stm_t_eval = stm_t(:,:,start_idx:target_idx);
% 
%     tcm_time = t(tcm_idx)';
% 
%     if start_idx == 1
%         vel_disp_flag = 0;
%         
%     else
%         vel_disp_flag = 1;
%     end

    

%     [~, tcm_dv_total] = calc_covariance_tcmdv(x(:), t_eval, t_s_eval, stm_t_eval, tcm_time, vel_disp_flag, P, simparams);
% 
% 
% 
%     tcm_time = traj.t(tcm_idx)';
% 
%     traj_eval.t = traj.t(start_idx:target_idx);
%     traj_eval.t_s = traj.t_s(start_idx:target_idx);
%     traj_eval.x_t = traj.x_t(start_idx:target_idx,:);
% 
% 
%     stm0M = invert_stm(traj.stm_t(:,:,start_idx), simparams);
%     traj_eval.stm_t = tmult(traj.stm_t(:,:,start_idx:target_idx), stm0M);
% 
% 
%     traj_eval.stm_t_i = traj.stm_t_i(start_node:target_node-1);
% 
% 
%     Qbarstart_tstart_t0 = traj.Q_t(:,:,start_node);
%     Qbar_t_tstart = tmult(traj_eval.stm_t, tmult(Qbarstart_tstart_t0, traj_eval.stm_t,[0 1]));
%     traj_eval.Q_t = traj.Q_t(:,:,start_idx:target_idx) - Qbar_t_tstart;
% 
%     traj_eval.Q_t_i = traj.Q_t_i(start_node:target_node-1);
% 
% 
%     maneuver_include = simparams.maneuverSegments>=traj.t_s(start_idx) & simparams.maneuverSegments < target_node;
%     deltaVs_nom_eval = deltaVs_nom(:,maneuver_include);
% 
% 
%     [~, tcm_dv_total] = calc_covariance_wQ_tcmdv(x, traj_eval, tcm_time, vel_disp_flag, deltaVs_nom_eval, P_i, simparams);

end