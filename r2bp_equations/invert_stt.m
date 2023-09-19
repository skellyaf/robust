function [istt] = invert_stt(stm_f, stt_f, simparams)
%invert_stt Takes a 2nd order state transition tensor in one direction (for
%example forward in time) and reverses it to go the opposite direction
%(like backwards in time)
%   Cite: Park, 2007, Einstein notation derivation (Eqn 15)
%   Inputs: 
%           stm_f: the forward stm for the corresponding timeframe
%           stt_f: the forward stt for the corresponding timeframe




% Calcuate the inverse stm for the corresponding timeframe
stm_r = invert_stm(stm_f,simparams);


% ATA = zeros(6,6,6);


T = stt_f;
A1 = stm_r;
A2 = stm_r;
% N=length(ATA);



% for i = 1:N
%     for j = 1:N
%         for k = 1:N
%             runsum=0;
%             for p = 1:N
%                 for q = 1:N
%                     runsum = runsum + T(i,p,q)*A1(p,j)*A2(q,k);
% %                     ATA(i,j,k) = ATA(i,j,k) + T(i,p,q)*A(p,j)*A(q,k);
%                 end
%             end
%             ATA(i,j,k) = runsum;
%         end
%     end
% end


%%%%%%%%%% VECTORIZED EQUIVALENT %%%%%%%%%%



pT = ipermute(T, [1, 4, 5, 2, 3]); %T
pA1 = ipermute(A1, [4, 2, 1, 3, 5]); %A1
pA2 = ipermute(A2, [5, 3, 1, 2, 4]); %A2

ATA = sum(  permute(  pT.*pA1.*pA2, [1, 2, 3, 4, 5]  ), [4, 5]   );

istt = -tmult(stm_r,ATA);



end