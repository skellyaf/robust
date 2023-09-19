function [traj_fd_new, tcm_idx_new] = replace_seg_fd(traj, traj_fd, tcm_idx, seg)
%replace_seg_fd Is used for finite difference testing/verification of
%analytical gradients. Specifically, when a segment duration is modified,
%all subsequent segments are starting at a different time, resulting in
%sensitivities that shouldn't exist. This function takes the trajectory
%parameters pre modification and only propagates/modifies the trajectory
%parameters for the current segment.

% List of trajectory parameters that need to be modified:

ppp=1;

% t
% t_s
% x_t
% x_i_f
% stm_t
% stm_t_i
% stm_i
% stt_t_i
% stt_i
% Q_i
% Q_t_i
% dQ_i
% dQ_t_i
% Q_t


% tcm_time(s) after seg?
% tcm_idx(s) after seg?

end