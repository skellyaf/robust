

A = randn(6,6);
B = randn(6,6,6);
C = randn(6,6);

testAC = zeros(6,6);

for i = 1:6
    for j  = 1:6
        for k = 1:6
            testAC(i,k) = testAC(i,k) + A(i,j)*C(j,k);
        end
    end
end


A*C - testAC
















T = randn(6,6,6);
A1 = randn(6,6);
A2 = A1;
N=6;


% 
% T = sym('T',[6 6 6], 'real');
% A1 = sym('A',[6 6], 'real');
% A2 = A1;

% e = ones(N,1);
% 
% 
% %%%%%%%%%% VECTORIZE THIS SOMETIME/SOMEHOW %%%%%%%%%%
% 
% for i = 1:N
%     for j = 1:N
%         for k = 1:N
%             runsum=0;
%             for p = 1:N
%                 for q = 1:N
% %                     runsum = runsum + A1(p,j)*A2(q,k);
%                     runsum = runsum + T(i,p,q)*A1(p,j)*A2(q,k);
%                 end
%             end
%             ATA(i,j,k) = runsum;
%             % the A part of runsum to this point
%             AAp = A1(:,j)*A1(:,k)';
% %             Apart = sum(sum(A1(:,j)*A1(:,k)'));
% %             simplify(runsum - Apart)
%             % Corresponds to the entirety of runsum:
%             M = squeeze(T(i,:,:)) .* AAp;
%             
%             ATA2(i,j,k) = sum(sum(squeeze(T(i,:,:)) .* AAp));
%             ATA3(i,j,k) = e' * M * e;
%         end
%     end
% end
% 
% 
% 
% 
% Tnew = zeros(6,6,6);
% 
% for i = 1:6
%     Tnew(:,:,i) = squeeze(T(i,:,:));
% end
% 
% 
% % tst = tmult(Tnew, A1*A1');
% 
% 
% 
% addpath('C:\Users\skell\Downloads\einsum.m')
% tsta = einsum(A1,A1,1,1);
% tstp = einsum(T,tsta,1,1)
% % tst = einsum(tstp,A1,3,1)



%%%% setup using ipermut and permute
%%% order: i (1), j (2), k (3), p (4), q (5)
% T( i(1),p(4),q(5) ) ---> 1, 4, 5, 2, 3
% A1( p(4), j(2) )
% A2 ( q(5), k(3) )

pT = ipermute(T, [1, 4, 5, 2, 3]);
pA1 = ipermute(A1, [4, 2, 1, 3, 5]);
pA2 = ipermute(A2, [5, 3, 1, 2, 4]);

pATA = sum(  permute(  pT.*pA1.*pA2, [1, 2, 3, 4, 5]  ), [4, 5]   );