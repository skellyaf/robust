function [Qbarf, dQbarf] = Qcombine(traj, idx_i, idx_f, simparams)
%Qcombine returns the cumulative Qbar from idx_i (assumes the initial
%Qbar is zero) to idx_f (initial to final)

% Notation: 
%%% Qbar is the accumulation of process noise covariance
%%% Qbarx (adding an identifier x) is the process noise covariance
%%% at location/time x
%%% Qbarx_t2_t1 is the contribution to Qbarx from t1 to t2 (backwards like
%%% a STM)

% Method: accumulate Qbarf as the sum of each of the involved segments'
% contributions. Compute in reverse so the STM to idx_f can be constructed
% sequentially. Like stm_f_j50; then stm_f_j50 * stm_j4f_j40; etc.

% Prep:
%%% 1. Find the initial and final segment numbers (segi and segf).
% Sub-steps:
%%% 1. Incorporate the Qbarf effect from the final segment. This is a
%%% direct extraction from Q_t_i{segment j} at the segment index for idx_f.
%%% 2. Incorporate the Qbarf effect from complete intermediate segments. Us
%%% the reverse sequentially computed stm to idx_f to then propagate Qbarj
%%% to Qbarf.
%%% 3. Check if the first segment (j is the segment index) is a fractional 
%%% segment. If so, computing the Qbarf contribution requires multiple
%%% steps as the propagations for Qbar all happened from the beginning of
%%% the segment (in Q_t_i cell structure).
%%%%% QbarjF_jf_i = QbarjF_jf_j0 - stm_jf_i * QbarjF_i_j0 * stm_jf_i'
%%% The fractional Qbar from segment j's effect on the final 

% simparams.dynSys = dynSys;


% Initial and final involved segment numbers
segi = traj.t_s(idx_i);
segf = traj.t_s(idx_f);

% Get the segment specific index for the first segment
segi_idx_i = find_cell_t_i_idx_2(traj.t,traj.t_s,idx_i, traj.stm_t_i, segi);

% If it is the last index in the specific segment, start at the beginning
% of the next segment.
fullInitialSeg = 0;
if segi_idx_i == size(traj.stm_t_i{segi},3) && segf > segi
    segi = segi + 1;
    fullInitialSeg = 1;
end




% Initialize structure for Qbarf
Qbarf = zeros(6,6);

if nargout == 2
    dQbarf = zeros(6,6,6);
end

iter = 1;
j = segf; % Start at the final segment and go backwards

while iter


    if j == segi && j == segf
        % Get the effect from the beginning of segi/segf (segj) to idx_f
        segj_idx_f = find_cell_t_i_idx_2(traj.t,traj.t_s,idx_f, traj.stm_t_i, j);
%         find_cell_t_i_idx(traj.t, traj.t_s, traj.stt_t_i, idx_f)

        if idx_f == 258
            ppp=1;
        end

        Qbarf_f_segj0 = traj.Q_t_i{j}(:,:,segj_idx_f);

        % Adding on dQ - mirroring Qbarf calculations
        if nargout == 2
            dQbarf_f_segj0 = traj.dQ_t_i{j}(:,:,:,segj_idx_f);
        end

        % If don't start from the beginning of the segment, subtract the
        % effect of Qbar on f from the beginning of the segment to
        % segf_idx_i

        segj_idx_i = find_cell_t_i_idx_2(traj.t,traj.t_s,idx_i, traj.stm_t_i, j);

        if fullInitialSeg && segj_idx_i == size(traj.Q_t_i{segi-1},3)
            segj_idx_i = 1;
        end


        if segj_idx_i > 1
            % The Qbar effect at i from the beginning of seg j to i
            Qbari_i_segj0 = traj.Q_t_i{segf}(:,:,segj_idx_i);
            % STMs required to propagate the above effect to f
            stm_i_j0 = traj.stm_t_i{j}(:,:,segj_idx_i);
            stm_j0_i = invert_stm(stm_i_j0, simparams);
            stm_f_j0 = traj.stm_t_i{j}(:,:,segj_idx_f);
            stm_f_i = stm_f_j0 * stm_j0_i;
            % The effect of Qbar at f from the beginning of segj to i
            Qbarf_i_segj0 = stm_f_i * Qbari_i_segj0 * stm_f_i';

            % The total Qbar effect at f: subtracting the effect from j0 to
            % i at f ... from the total effect from j0 to f at f
            Qbarf = Qbarf_f_segj0 - Qbarf_i_segj0;

            % dQbar calculations
            if nargout == 2
                stt_f_j0 = traj.stt_t_i{j}(:,:,:,segj_idx_f);
                stt_i_j0 = traj.stt_t_i{j}(:,:,:,segj_idx_i);
                stt_j0_i = invert_stt(stm_i_j0, stt_i_j0, simparams);

                stt_f_i = tensorCombine(stm_j0_i, stt_j0_i, stm_f_j0, stt_f_j0);

                dQbar_i_segj0 = traj.dQ_t_i{j}(:,:,:,segj_idx_i);

                dQbarf = dQseparate_prior(dQbarf_f_segj0, dQbar_i_segj0, Qbari_i_segj0, stm_i_j0, stm_f_i, stt_f_i, simparams);
            end
        else
            Qbarf = Qbarf_f_segj0;

            % dQbar
            if nargout == 2
                dQbarf = dQbarf_f_segj0;
            end
        end

        iter = 0;
    elseif j == segf
        % Include the effect from the beginning of segf to idx_f
        segf_idx_f = find_cell_t_i_idx_2(traj.t,traj.t_s,idx_f, traj.stm_t_i, segf);

        Qbarf_f_segf0 = traj.Q_t_i{segf}(:,:,segf_idx_f);

        Qbarf = Qbarf + Qbarf_f_segf0;

        % dQbar
        if nargout == 2
            dQbarf = traj.dQ_t_i{segf}(:,:,:,segf_idx_f);
            % Also get the reverse STT
            stt_f_Cidx = traj.stt_t_i{segf}(:,:,:,segf_idx_f);
        end

        % Get the STM from the beginning of seg f to idx_f. This starts the
        % STM reverse combination. Cidx stands for current index.
        stm_f_Cidx = traj.stm_t_i{segf}(:,:,segf_idx_f);

       
    elseif j > segi
        % Include the entire effect of the full segment between segf & segi
        Qbarj = traj.Q_i(:,:,j);
        
        Qbarf_jf_j0 = stm_f_Cidx * Qbarj * stm_f_Cidx';

        Qbarf = Qbarf + Qbarf_jf_j0;

        stm_jf_j0 = traj.stm_t_i{j}(:,:,end); % STM from the beginning to end of the current segment, j

        % Include the effect of dQ also
        if nargout == 2
            % Assemble the pieces to combine dQ_jf_j0 with dQbarf
            dQbar_jf_j0 = traj.dQ_t_i{j}(:,:,:,end);
            dQbarf = dQcombine_sequential(dQbar_jf_j0, dQbarf, Qbarj, stm_jf_j0, stm_f_Cidx, stt_f_Cidx);
        end


        
        
        % Assemble the reverse total STT
        if nargout == 2
            stt_jf_j0 = traj.stt_t_i{j}(:,:,:,end);
            stt_f_Cidx = tensorCombine(stm_jf_j0, stt_jf_j0, stm_f_Cidx, stt_f_Cidx);
        end

        % Continue assembling the reverse total STM
        stm_f_Cidx = stm_f_Cidx * stm_jf_j0; % STT from the beginning to end of the current segment, j 



    elseif j == segi
        % Check if the whole segment or a partial segment contribution is
        % being included

        if fullInitialSeg
            % Include the entire effect of the full segment between segf &
            % segi (same as above full segment)
            Qbarj = traj.Q_i(:,:,j);
            
            Qbarf_jf_j0 = stm_f_Cidx * Qbarj * stm_f_Cidx';
    
            Qbarf = Qbarf + Qbarf_jf_j0;

            % Include the effect of dQ also
            if nargout == 2
                stm_jf_j0 = traj.stm_t_i{j}(:,:,end); % STM from the beginning to end of the current segment, j
                % Assemble the pieces to combine dQ_jf_j0 with dQbarf
                dQbar_jf_j0 = traj.dQ_t_i{j}(:,:,:,end);
                dQbarf = dQcombine_sequential(dQbar_jf_j0, dQbarf, Qbarj, stm_jf_j0, stm_f_Cidx, stt_f_Cidx);
            end
            
        else
            % Only include a fraction of the initial segment
%             QbarjF_jf_i = QbarjF_jf_j0 - stm_jf_i * QbarjF_i_j0 * stm_jf_i'
            % Get the segment index corresponding to idx_i
            segj_idx_i = find_cell_t_i_idx_2(traj.t,traj.t_s,idx_i, traj.stm_t_i, j);
            % Get the STM from idx_i to the end of the segment
            stm_i_j0 = traj.stm_t_i{j}(:,:,segj_idx_i);
            stm_j0_i = invert_stm(stm_i_j0, simparams);
%             stm_jf_j0 = traj.stm_i(:,:,j);
            stm_jf_j0 = traj.stm_t_i{j}(:,:,end);
            stm_jf_i = stm_jf_j0 * stm_j0_i;

            % Get Qbar at the end of segment j (jf) for the entire segment
            Qbarjf_jf_j0 = traj.Q_i(:,:,j);

            % Get the Qbar at i from the effect j0 to i
            Qbari_i_j0 = traj.Q_t_i{j}(:,:,segj_idx_i);
            % Propagate that to the end of segment j (jf)
            Qbarjf_i_j0 = stm_jf_i * Qbari_i_j0 * stm_jf_i';

            % Subtract the impact from j0 to i at jf from the entire
            % segment's impact at jf
            Qbarjf_jf_i = Qbarjf_jf_j0 - Qbarjf_i_j0;

            % Propagate that to idx_f
            Qbarf_jf_i = stm_f_Cidx * Qbarjf_jf_i * stm_f_Cidx';

            % Add it to Qbarf
            Qbarf = Qbarf + Qbarf_jf_i;

            % Do the same with dQ
            if nargout == 2
                % Remove dQ from the beginning of the segment to the
                % starting index from dQ for the entire segment
                dQbar_jf_j0 = traj.dQ_t_i{j}(:,:,:,end); % dQbar for the entire segment 
                dQbar_i_j0 = traj.dQ_t_i{j}(:,:,:,segj_idx_i); % dQbar from the beginning of the segment to the start idx

                % Assemble STT from i to jf
                stt_jf_j0 = traj.stt_t_i{j}(:,:,:,end);
                stt_i_j0 = traj.stt_t_i{j}(:,:,:,segj_idx_i);
                stt_j0_i = invert_stt(stm_i_j0, stt_i_j0, simparams);
                stt_jf_i = tensorCombine(stm_j0_i, stt_j0_i, stm_jf_j0, stt_jf_j0);
                dQbar_jf_i = dQseparate_prior(dQbar_jf_j0, dQbar_i_j0, Qbari_i_j0, stm_i_j0, stm_jf_i, stt_jf_i, simparams);

                % Combine dQbar_jf_i with all of the assembled pieces after it until f
                dQbarf = dQcombine_sequential(dQbar_jf_i, dQbarf, Qbarjf_jf_i, stm_jf_i, stm_f_Cidx, stt_f_Cidx);
            end
            
            ppp=1;
        end



        % End loop if at the initial segment
        iter = 0;
    else
        assert(0,'Something went wrong.');
    end




    

    % Decrement j
    j = j - 1;
end


end