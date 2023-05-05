function tcm_dv_total = calc_tcm_dv(tcm_idx, x, t, t_s, stm_t, P, start_idx, target_idx, simparams)

    t_eval = t(start_idx:target_idx);
    t_s_eval = t_s(start_idx:target_idx);
    stm_t_eval = stm_t(:,:,start_idx:target_idx);

    tcm_time = t(tcm_idx)';

    if start_idx == 1
        vel_disp_flag = 0;
        
    else
        vel_disp_flag = 1;
    end

    

    [~, tcm_dv_total] = calc_covariance_tcmdv(x(:), t_eval, t_s_eval, stm_t_eval, tcm_time, vel_disp_flag, P, simparams);

end