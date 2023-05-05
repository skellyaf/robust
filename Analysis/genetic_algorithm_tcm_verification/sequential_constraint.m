function [cin, ceq] = sequential_constraint(tcm_idx, x, t, t_s, stm_t, P, start_idx, target_idx, simparams)

    ceq = [];
%     cin = zeros(1,length(tcm_idx) - 1);
%     
%     for i = 1:length(tcm_idx) - 1
%         cin(i) = tcm_idx(i) - tcm_idx(i+1) + 1; % plus one bc it is less than or equal to, and they're integers
%     end




    t_eval = t(start_idx:target_idx);
    t_s_eval = t_s(start_idx:target_idx);
    stm_t_eval = stm_t(:,:,start_idx:target_idx);

    tcm_time = t(tcm_idx)';

    if start_idx == 1
        vel_disp_flag = 0;
        
    else
        vel_disp_flag = 1;
    end





    [P, tcm_dv_total] = calc_covariance_tcmdv(x(:), t_eval, t_s_eval, stm_t_eval, tcm_time, vel_disp_flag, P, simparams);

    cin = sqrt(trace(P(1:3,1:3))) - simparams.P_max_r;

end