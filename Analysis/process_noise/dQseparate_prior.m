function [dQ_2_1] = dQseparate_prior(dQ_2_0, dQ_1_0, Q_1_0, stm_1_0, stm_2_1, stt_2_1, simparams)
%dQseparate_prior removes the dQ effect at the beginning of a segment and
%returns the custom dQ portion toward the end of a segment. 
%Ex: Have dQ2_2_0 and dQ1_1_0 and want dQ2_2_1

stm_0_1 = invert_stm(stm_1_0, simparams);

dphi21_Q10_phi21_dx0 = einsum3(dQ_1_0, stm_2_1, stm_2_1, [4, 5, 3], [1, 4], [2, 5], 3);

dQ2_2_0_dx1 = einsum2(dQ_2_0 - dphi21_Q10_phi21_dx0, stm_0_1, [1, 2, 4], [4, 3], 3);

t1 = einsum3(Q_1_0, stt_2_1, stm_2_1, [4, 5], [1, 4, 3], [2, 5], 3);
t2 = einsum3(Q_1_0, stm_2_1, stt_2_1, [4, 5], [1, 4], [2, 5, 3], 3);

dQ_2_1 = dQ2_2_0_dx1 - t1 - t2;

end