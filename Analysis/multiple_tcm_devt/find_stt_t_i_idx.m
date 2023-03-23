function [stt_idx, seg] = find_stt_t_i_idx(t, t_s, stt_t_i, t_idx)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
seg = t_s(t_idx);
cum_pre_idx = 0;
if seg > 1
    for i = 1:seg - 1
        if i == 1
            nAdd = 0;
        else
            nAdd = 1;
        end
    
        cum_pre_idx = cum_pre_idx + size(stt_t_i{i},4) - nAdd;
    
    
    end
    stt_idx = t_idx - cum_pre_idx + 1;
else
    stt_idx = t_idx;



end