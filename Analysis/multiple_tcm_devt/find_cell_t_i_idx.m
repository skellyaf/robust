function [stt_idx, seg] = find_cell_t_i_idx(t, t_s, stt_t_i, t_idx)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%   CAN ACCEPT STM_T_I OR STT_T_I STRUCTURE


seg = t_s(t_idx);
cum_pre_idx = 0;
i=1;
dims=1;
if seg > 1
    for i = 1:seg - 1
        if i == 1
            nAdd = 0;
        else
            nAdd = 1;
        end
    
        dims = size(stt_t_i{i});
        cum_pre_idx = cum_pre_idx + dims(end) - nAdd;
    
    
    end
    stt_idx = t_idx - cum_pre_idx + 1;
else
    stt_idx = t_idx;



end