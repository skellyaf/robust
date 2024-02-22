function [dstmCkClastdxi, dstmCkClastddti] = calc_dstmCkClast(x, x_i_f, t, t_s, stm_t, stm_i, stm_t_i, stt_t_i, tCk, tClast, i, simparams)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% Symplectic unit matrix
% J = [zeros(3,3), eye(3,3); -eye(3,3), zeros(3,3)];
mu = simparams.mu;
m = simparams.m;
n = simparams.n;
nsv = simparams.nsv;

% The initial time for the current segment
if i == 1
    t0_i = 0;
else
    t0_i = sum(x(m,1:i-1));
end

tf_i = sum(x(m,1:i)); % The final time for the current segment


tCk_idx = find(t==tCk);
tClast_idx = find(t==tClast);


%%%%% 6 LOGIC CASES FOR dstmCkClastdxi

%% Logic case 1: [t0,i < tf,i < tClast < tCk]
if tf_i <= tClast
    % State adjustments happen before any of the segments involved in the corrections 
    dstmCkClastdxi = zeros(nsv,nsv,nsv);
    dstmCkClastddti = zeros(nsv,nsv);

%% Logic case 2: [t0,i < tClast < tf,i < tCk] (TESTED - IS WORKING!)
elseif t0_i <= tClast && tf_i <= tCk

    ip_idx = find(t==tf_i);
    stmCkip = dynCellCombine(t, t_s, ip_idx, tCk_idx, simparams, stm_t_i, stt_t_i);

    i0_idx = find(t==t0_i);
    [stmClasti0, sttClasti0] = dynCellCombine(t, t_s, i0_idx, tClast_idx, simparams, stm_t_i, stt_t_i);
%     stmi0Clast = -J * stmClasti0' * J;
    stmi0Clast = invert_stm(stmClasti0, simparams);

    % Using matrix inverse property
    dstmClasti0_inv = -tmult(stmi0Clast,  tmult(sttClasti0, stmi0Clast)  );

    dstmCkClastdxi = tmult( stmCkip, tmult(  stm_i(:,:,i), dstmClasti0_inv  ) + tmult(  stt_t_i{i}(:,:,:,end), stmi0Clast  )   );

    % Calc dstmCkClastddti
    stmifClast = dynCellCombine(t, t_s, tClast_idx, ip_idx, simparams, stm_t_i, stt_t_i);
    
%     dstmCkClastddti = stmCkip * r2bp_A_matrix(x_i_f(:,i),mu) * stmifClast;
    dstmCkClastddti = stmCkip * A_matrix(x_i_f(:,i),simparams) * stmifClast;
    

%% Logic case 3: [t0,i < tClast < tCk < tf,i] (TESTED - IS WORKING!)
elseif t0_i <= tClast && tCk < tf_i

    i0_idx = find(t==t0_i);
    [stmCki0, sttCki0] = dynCellCombine(t, t_s, i0_idx, tCk_idx, simparams, stm_t_i, stt_t_i);  
    
    [stmClasti0, sttClasti0] = dynCellCombine(t, t_s, i0_idx, tClast_idx, simparams, stm_t_i, stt_t_i);
%     stmi0Clast = -J * stmClasti0' * J;
    stmi0Clast = invert_stm(stmClasti0, simparams);

    % Using matrix inverse property
    dstmClasti0_inv = -tmult(stmi0Clast,  tmult(sttClasti0, stmi0Clast)  );

    dstmCkClastdxi = tmult( stmCki0, dstmClasti0_inv ) + tmult( sttCki0, stmi0Clast );

    % Calc dstmCkClastddti
    dstmCkClastddti = zeros(nsv,nsv);



%% Logic case 4: [tClast < t0,i < tCk < tf,i] (TESTED - IS WORKING!)
elseif t0_i <= tCk && tCk < tf_i



    i0_idx = find(t==t0_i);

    [~, sttCki0] = dynCellCombine(t, t_s, i0_idx, tCk_idx, simparams, stm_t_i, stt_t_i);

    stmi0Clast = dynCellCombine(t, t_s, tClast_idx, i0_idx, simparams, stm_t_i);

    dstmCkClastdxi = tmult( sttCki0, stmi0Clast );

    % Calc dstmCkClastddti
    dstmCkClastddti = zeros(nsv,nsv);


%% Logic case 5: [tClast < tCk < t0,i < tf,i]
elseif tCk < t0_i
    % State adjustments happen after any of the segments involved in the corrections 
    dstmCkClastdxi = zeros(nsv,nsv,nsv);

    % Calc dstmCkClastddti
    dstmCkClastddti = zeros(nsv,nsv);

%% Logic case 6: [tClast < t0,i < tf,i < tCk]
elseif tf_i <= tCk
    i0_idx = find(t==t0_i);
    ip_idx = find(t==tf_i);

    stmCkip = dynCellCombine(t, t_s, ip_idx, tCk_idx, simparams, stm_t_i);

    stmi0Clast = dynCellCombine(t, t_s, tClast_idx, i0_idx, simparams, stm_t_i);

    stti = stt_t_i{i}(:,:,:,end);
    stmi = stm_t_i{i}(:,:,end);

    dstmCkClastdxi = tmult(stmCkip, tmult(stti, stmi0Clast));

    dstmCkClastddti = stmCkip * A_matrix(x_i_f(:,i),simparams) * stmi * stmi0Clast;

%% Something is messed up if we're here
else
    assert(0,'None of the logic cases were satisfied, something is amiss!');

end





end