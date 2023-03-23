function xstmdot = r2bp_stt2_de(t, x, mu)

% global mu;

r = [x(1); x(2); x(3)];
v = [x(4); x(5); x(6)];
stm = reshape(x(7:42), 6,6);

% ADDED STT
stt2 = reshape(x(43:258),6,6,6);


f = [zeros(3,3), eye(3,3); 
    -mu/norm(r)^3 * eye(3,3), zeros(3,3)];

xdot = f * [r; v];

% Partial wrt state vector
A = r2bp_A_matrix([r; v], mu);

% Second partial wrt state vector (returns a tensor 6x6x6)
Fab = r2bp_stt2_tensor([r;v], mu);

% STM propagation equation
stmdot = A*stm;

% STT propagation equation
% stt2dot = zeros(6,6,6);
% for i = 1:6
%     stt2dot(:,:,i) = A*stt2(:,:,i) + Fab(:,:,i)*stm*stm;
% end



stt2dot = tensorCombine(stm,stt2,A,Fab);

% Return
xstmdot = [xdot; stmdot(:); stt2dot(:)];
        

end





% 
% 
% A = rand(4,10,3);
% B = rand(10,16);
% 
% C(:,:,1) = A(:,:,1)*B;
% C(:,:,2) = A(:,:,2)*B;
% C(:,:,3) = A(:,:,3)*B;
% >> [m,n,r] = size(A);
% out = permute(reshape(reshape(permute(A,[1 3 2]),[],n)*B,m,r,[]),[1 3 2]);