function [dPCkminusd] = calc_dPckMinus(stmCkClast, dstmCkClast, R, Tlast, dTlast, PClast_minus, dPClast_minus)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here




% % % % if k == 1
% % % % %     PClast_minus = simparams.P_initial;
% % % % 
% % % %     dPCkminusd = tmult(stmCkClast*PClast_minus,dstmCkClast,[0 1])         + tmult(dstmCkClast,PClast_minus*stmCkClast');
% % % % 
% % % % else

%     R = simparams.R;

    G = [zeros(3,3); eye(3,3)]; % mapping matrix to vv submatrix

    
    INlast = eye(6) + [zeros(3,6); Tlast];
    
    depth = size(dTlast,3);

    dINlast = zeros(6,6,depth);
    dINlast(4:6,:,:) = dTlast;
    
    
    t1 = tmult(stmCkClast*G*R*G',dstmCkClast,[0 1])                            + tmult(dstmCkClast,G*R*G'*stmCkClast');
    
    t2 = tmult(stmCkClast*INlast*PClast_minus*INlast',dstmCkClast,[0 1])         + tmult(dstmCkClast,INlast*PClast_minus*INlast'*stmCkClast');
    
    t3 = tmult(stmCkClast,tmult(dINlast,PClast_minus*INlast'*stmCkClast'))          + tmult(stmCkClast*INlast*PClast_minus,tmult(dINlast,stmCkClast,[1 1]));
    
    t4 = tmult(stmCkClast*INlast,tmult(dPClast_minus,INlast'*stmCkClast'));
    
    
    dPCkminusd = t1 + t2 + t3 + t4;

% % % % end


end