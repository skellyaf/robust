function dT = T_partial(stmNI, STTnCdx)
% The partial derivative of a corrective maneuver with respect to either
% states or times
% Really just vector, matrix, and tensor math to make other functions
% easier to implement
% First order linear targeting


depth = size(STTnCdx,3);
nsv = size(stmNI,2);
T_analytical_partial = zeros(3,nsv,depth);



%%%%%% NOT SURE THIS IS THE RIGHT WAY TO DEAL WITH THIS
if rank(stmNI(1:3,4:6)) > 0
    invstmNIrv = inv(stmNI(1:3,4:6));


        lhs = tmult( invstmNIrv ,  tmult(   STTnCdx(1:3,4:6,:)    , invstmNIrv * stmNI(1:3,1:3),[0 0]) ,[0 0]) - tmult(invstmNIrv, STTnCdx(1:3,1:3,:),[0 0]);
        T_analytical_partial(1:3,1:3,:) = lhs;
end

dT = T_analytical_partial;


end