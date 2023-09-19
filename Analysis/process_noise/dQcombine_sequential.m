function [dQ_2_0] = dQcombine_sequential(dQ_1_0, dQ_2_1, Q_1_0, stm_1_0, stm_2_1, stt_2_1)
%dQcombine_sequential combines two sequential process noise covariance
%state tensors into one that grew from the beginning to the end of the time
%period. 



dphi21_Q10_phi21_dx0 = einsum3(dQ_1_0, stm_2_1, stm_2_1, [4, 5, 3], [1, 4], [2, 5], 3);


t1 = einsum3(Q_1_0, stt_2_1, stm_2_1, [4, 5], [1, 4, 3], [2, 5], 3);
t2 = einsum3(Q_1_0, stm_2_1, stt_2_1, [4, 5], [1, 4], [2, 5, 3], 3);

dQ2_2_1_dx0 = einsum2(dQ_2_1 + t1 + t2, stm_1_0, [1, 2, 4], [4, 3], 3);


dQ_2_0 = dQ2_2_1_dx0 + dphi21_Q10_phi21_dx0;


end