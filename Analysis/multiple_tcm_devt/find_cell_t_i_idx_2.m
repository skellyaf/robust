function [seg_index, seg] = find_cell_t_i_idx_2(t, t_s, t_idx, stm_t_i, desired_seg)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



seg = t_s(t_idx);
idxs = find(seg==t_s);
seg_index = find(idxs == t_idx);

if seg > 1
    seg_index = seg_index + 1;
end



if seg ~= desired_seg 
    if seg_index == size(stm_t_i{seg},3)
        seg = seg+1;
        seg_index = 1;
    else
        assert(0,['Error!'])
    end
end






% seg = t_s(t_idx);
% cum_pre_idx = 0;
% i=1;
% dims=1;
% if seg > 1
%     for i = 1:seg - 1
%         if i == 1
%             nAdd = 0;
%         else
%             nAdd = 1;
%         end
%     
%         dims = size(stt_t_i{i});
%         cum_pre_idx = cum_pre_idx + dims(end) - nAdd;
%     
%     
%     end
%     stt_idx = t_idx - cum_pre_idx + 1;
% else
%     stt_idx = t_idx;



end