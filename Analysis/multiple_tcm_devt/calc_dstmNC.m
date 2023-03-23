function [dstmNCdxi,dstmNCddti] = calc_dstmNC(x, x_i_f, t, t_s, stm_t, stm_i, stm_t_i, stt_t_i, tCk, tCk_idx, tCk_seg, stmN0, i, simparams)

mu = simparams.mu;

% Symplectic unit matrix
% J = [zeros(3,3), eye(3,3); -eye(3,3), zeros(3,3)];

%% 4 LOGIC CASES 
%%
%% Logic case 1: segment i is before the segment containing the TCM Ck 
if i < tCk_seg

    % ---- dstmNCdxi
    dstmNCdxi = zeros(6,6,6);

    % ---- dstmNCddti
    dstmNCddti = zeros(6,6);


%% Logic case 2: segment i contains the TCM
elseif i == tCk_seg

    % Get the time the correction segment starts = t_corrSegStart
    if tCk_seg > 1
        t_corrSegStart = sum(x(7,1:tCk_seg-1));
    else
        t_corrSegStart = 0;
    end
    
    % Get the time the correction segment ends = t_corrSegEnd
    t_corrSegEnd = sum(x(7,1:tCk_seg));

    % Get stmNip0
    ip0_idx = find(t==t_corrSegEnd);
    target_idx = length(t);
    stmNip0 = dynCellCombine(t, t_s, ip0_idx, target_idx, simparams, stm_t_i);

    % Get sttCi0 and stmi0C
    i0_idx = find(t==t_corrSegStart);
    [stmCi0, sttCi0] = dynCellCombine(t, t_s, i0_idx, tCk_idx, simparams, stm_t_i, stt_t_i);
%     stmi0C = -J * stmCi0' * J;
    stmi0C = invert_stm(stmCi0, simparams);
    % Using matrix inverse property
    dstmCi0_inv = -tmult(stmi0C,  tmult(sttCi0, stmi0C)  );

    % calc dstmNCdxi
    dstmNCdxi = tmult(stmNip0,   tmult(stt_t_i{i}(:,:,:,end), stmi0C)   + tmult(stm_i(:,:,i), dstmCi0_inv)       );

    % Get stmifC
    stmifC = dynCellCombine(t, t_s, tCk_idx, ip0_idx, simparams, stm_t_i);

    % ---- dstmNCddti
    % Need state at end of correction segment
    x_c_f = x_i_f(:,i);
%     dstmNCddti = stmNip0 * r2bp_A_matrix(x_c_f, mu) * stmifC;
    dstmNCddti = stmNip0 * stmDot(x_c_f, simparams) * stmifC;

%% Logic case 3: segment i is after the segment containing the TCM
% but not after the target
elseif t(end) > sum(x(7,1:i-1))

    % Get the time elements for segment i
    t_i = t(t_s==i);

    % Get dstmNif
    target_idx = length(t);
    if_idx = find(t==t_i(end));

    stmNif = dynCellCombine(t, t_s, if_idx, target_idx, simparams, stm_t_i);

    % Get smi0C
    if i == 1
        i0_idx = 1;
    else
        i0_idx = find(t==t_i(1))-1;
    end

    stmi0C = dynCellCombine(t, t_s, tCk_idx, i0_idx, simparams, stm_t_i);

    % ---- dstmNCdxi
    dstmNCdxi = tmult(stmNif, tmult(stt_t_i{i}(:,:,:,end), stmi0C));

    % Get stmifC
    stmifC = dynCellCombine(t, t_s, tCk_idx, if_idx, simparams, stm_t_i);

    % ---- dstmNCddti
%     dstmNCddti = stmNif * r2bp_A_matrix(x_i_f(:,i), mu) * stmifC;
    dstmNCddti = stmNif * stmDot(x_i_f(:,i), simparams) * stmifC;

%% Logic case 4: segment i is after the target (coast segment)
else 
    % ---- dstmNCdxi
    dstmNCdxi = zeros(6,6,6);
    % ---- dstmNCddti
    dstmNCddti = zeros(6,6);


end




end