function [stm_combined, stt_combined] = dynCellCombine(t, t_s, total_idx_start, total_idx_end, simparams, stm_t_i, stt_t_i)
%dynCellCombine combines multiple STMs and STTs in the cell structures stm_t_i and stt_t_i to
%create a single STT that spans from t_idx_start to t_idx_end
%(stt_end_start)

if nargout == 2 && nargin == 6
    assert(0,'You forgot to input the STT cell structure!');
end

if length(total_idx_start) > 1
    total_idx_start = total_idx_start(1);
end
if length(total_idx_end) > 1
    total_idx_end = total_idx_end(1);
end


if nargout == 2
    stt_flag = 1;
else
    stt_flag = 0;
end
nsv = simparams.nsv;
% Enable the abilty to request an inverse STM or STT by providing a start
% index that is greater than the end index
if total_idx_start > total_idx_end
    invert_at_end = 1;
    save_total_idx_end = total_idx_end;
    total_idx_end = total_idx_start;
    total_idx_start = save_total_idx_end;
    skip = 0;
elseif total_idx_start < total_idx_end
    invert_at_end = 0;
    skip = 0;
else
    skip = 1;
    stm_combined = eye(nsv);
    stt_combined = zeros(nsv,nsv,nsv);
end

% STM and STT indices indicate the end of a segment, but we want to start
% at the beginning of a segment to go from beginning to end time so we need
% to subtract one, unless it is already starting from the beginning (1)

if ~skip
    t_start = t(total_idx_start);
    t_end = t(total_idx_end);
    
    seg_t_end = t_s(total_idx_end);
    seg_t_start = t_s(total_idx_start);
    
    


    % If the requested index is the final index of a segment, just start
    % from the beginning of the next segment from the cell structure
    % stm_t_i and stt_t_i
    seg_t_start_incides = find(t_s == seg_t_start);
    if total_idx_start == seg_t_start_incides(end) && seg_t_end > seg_t_start
        seg_t_start = seg_t_start + 1;
    end

    
    
    for i = seg_t_start:seg_t_end
        seg_times = t(t_s==i);


        seg_idx = find(t_s==i);
            
        % Get seg_idx_start
        seg_idx_start = seg_idx(1);
        seg_idx_end = seg_idx(end);

%         if t(total_idx_end) < seg_times(end)
        if total_idx_end < seg_idx_end
            % Don't go all the way to the end of the current cell

            end_logical = t(t_s==i)==t(total_idx_end);
            
%             if length(end_logical) > 1
%                 end_logical = end_logical(1);
%             end


            %%%%% TODO: WHAT WAS THIS NEXT SECTION OF CODE FOR? THERE
            %%%%% APPEARS TO BE A BUG IN THE LOGIC OF THE SIZE COMPARISON
            %%%%% IN THE NEXT IF STATEMENT....COMMENTING OUT FOR NOW TO
            %%%%% SEE IF IT CAUSES AN ERROR SOMEWHERE ELSE.



            %%%%%%% POTENTIAL ANSWER: I THINK IT WAS FOR THE FACT THAT t_s
            %%%%%%% IS STORED WITH A DIFFERENT NUMBER OF INDICES THAN
            %%%%%%% stm_t_i{} CONTAINS. t and t_s start at 0, then go to 1,
            %%%%%%% 2, 3, 4 sequentially regardless of whether or not it is
            %%%%%%% a new trajectory segment or not. stm_t_i, on the other
            %%%%%%% hand, effectively starts at 0 at the beginning of each
            %%%%%%% segment.

            if size(end_logical,1) ~= size(stm_t_i{i},3)
                end_logical = [false; end_logical];
            end

%%%%%%%%%%%%%%%%%%%----------------------------




            if stt_flag
                stt_if_i0 = stt_t_i{i}(:,:,:,end_logical);
            end
            stm_if_i0 = stm_t_i{i}(:,:,end_logical);
        else
            % Go all the way to the end of the current cell
            if stt_flag
                stt_if_i0 = stt_t_i{i}(:,:,:,end);
            end
            stm_if_i0 = stm_t_i{i}(:,:,end);
        end
    



        
        % If the requested start index is after the beginning of the current segment
        if total_idx_start > seg_idx_start
%         if t(total_idx_start) > seg_times(1)

            



            % Then need to combine from the beginning of the segment to the end index      
            cell_idx_start = find_cell_t_i_idx(t,t_s,stm_t_i,total_idx_start);
    
            stm_i_i0 = stm_t_i{i}(:,:,cell_idx_start);
            stm_i0_i = invert_stm(stm_i_i0, simparams);
    
            stm_if_i = stm_if_i0 * stm_i0_i;
            
            if stt_flag
                % Then STTs
                stt_i_i0 = stt_t_i{i}(:,:,:,cell_idx_start);
                
                % Invert stt_i_0 to combine in subsequent step
                stt_i0_i = invert_stt(stm_i_i0,stt_i_i0,simparams);
        
                % Then combine STTs to get stt_if_i
                stt_if_i = tensorCombine(stm_i0_i, stt_i0_i, stm_if_i0, stt_if_i0);
            end
            
        else
            % If it starts from the beginning of the segment, only need to go to the end index
            stm_if_i = stm_if_i0;
            if stt_flag
                stt_if_i = stt_if_i0;
            end
    
        end
    
        % Now combine
        if i == seg_t_start
            % Initialize
            stm_combined = stm_if_i;
            if stt_flag
                stt_combined = stt_if_i;
            end
        else
            % Sequentially combine
            if stt_flag
                stt_combined = tensorCombine(stm_combined, stt_combined, stm_if_i, stt_if_i);
            end
            stm_combined = stm_if_i * stm_combined;
        end
    
    
    
    end
    
    
    if invert_at_end
        if stt_flag
            stt_combined = invert_stt(stm_combined, stt_combined, simparams);
        end
        stm_combined = invert_stm(stm_combined, simparams);
    
    end
end





end