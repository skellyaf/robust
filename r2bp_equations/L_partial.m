function dL = L_partial(W, T, dstmNC, dT)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dW = zeros(3,6,size(dstmNC,3));
dL = dW;
dstmNC_vr = dstmNC(4:6,1:3,:);
dstmNC_vv = dstmNC(4:6,4:6,:);
dW(:,1:3,:) = dstmNC_vv;
dW(:,4:6,:) = -dstmNC_vr;


dL(:,1:3,:) = tmult(dW,T,[0 1]) + tmult(W,dT,[0 1]);


end