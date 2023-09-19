function [stt12] = tensorCombine(stm1, stt1, stm2, stt2)


ATA = zeros(6,6,6);


% T = stt2;
% A1 = stm1;
% A2 = stm1;
% N=length(ATA);
% 
% 
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


%%%% setup using ipermute and permute from the help in this stack overflow thread: https://stackoverflow.com/questions/55913093/element-wise-matrix-multiplication-for-multi-dimensional-array 
%%% order: i (1), j (2), k (3), p (4), q (5)
% T( i(1),p(4),q(5) ) ---> 1, 4, 5, 2, 3
% A1( p(4), j(2) )
% A2 ( q(5), k(3) )
% 

pT = ipermute(stt2, [1, 4, 5, 2, 3]); %T
pA1 = ipermute(stm1, [4, 2, 1, 3, 5]); %A1
pA2 = ipermute(stm1, [5, 3, 1, 2, 4]); %A2

ATA = sum(  permute(  pT.*pA1.*pA2, [1, 2, 3, 4, 5]  ), [4, 5]   );
stt12 = ATA + tmult(stm2,stt1,[0 0]);


end