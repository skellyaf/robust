function [Qbar_events, Qbar_t] = Qcombine_2(x, traj, tcm_time, simparams)
%Qcombine returns Qbar at each increment along a nominal trajectory,
%considering all events/incorporations into the state dispersion.
%This 2nd version returns Qbar at each increment along the trajectory so
%only one set of Qbar calculations need to be performed for each nominal
%trajectory.

% Notation: 
%%% Qbar is the accumulation of process noise covariance
%%% Qbarx (adding an identifier x) is the process noise covariance
%%% at location/time x
%%% Qbarx_t2_t1 is the contribution to Qbarx from t1 to t2 (backwards like
%%% a STM)

% Method: accumulate Qbarf as the sum of each of the involved segments'
% contributions. 


%% Code start

% Get event info
[event_times, event_indicator] = define_events_v2(x(:), traj.t, tcm_time, simparams);



event_idx_logical = logical(sum(traj.t'==event_times', 1));    
event_idxs = find(event_idx_logical);

if nargout > 1
    Qbar_t = zeros(6,6,length(traj.t)); % for storing at each time index
end




Qbar_events = zeros(6,6,length(unique(event_times))); % if wanting to store only for each event, maybe a better idea

for i = 1:length(event_times)
    if i == 1
        % Starting from the beginning of the trajectory until the first
        % event, Qbart is exactly as propagated in Q_t

        
        Qbar_events(:,:,1) = traj.Q_t(:,:,event_idxs(1));



        if nargout > 1
            Qbar_t(:,:,1:event_idxs(1)) = traj.Q_t(:,:,1:event_idxs(1));
        end



    else

        % After the first event (after the Qbart propagation from the
        % beginning of the trajectory has been included in a dispersion
        % covariance correction calculation), Qbar includes Q_t and
        % subtracts the Q_t contribution prior to the previous event.
   

        stm_ei_0 = traj.stm_t(:,:,event_idxs(i));

        stm_eim_0 = traj.stm_t(:,:,event_idxs(i-1));
        stm_0_eim = invert_stm(stm_eim_0, simparams);
        stm_ei_eim = stm_ei_0 * stm_0_eim;

        Qbar_events(:,:,i) = traj.Q_t(:,:,event_idxs(i)) - stm_ei_eim * traj.Q_t(:,:,event_idxs(i-1)) * stm_ei_eim';



        if nargout > 1

            % Weird STM notation below. The first part, stm_t, means it is
            % a stm time history tensor, not just a single STM. The second
            % piece of data (stm_t_eim) indicates where the time index
            % begins. The third and fourth pieces of data (stm_t_eim_ei_0)
            % indicate the tranditional stm endpoints (reference points for
            % the dynamics/partial derivatives).

            stm_t_eim_ei_0 = traj.stm_t(:,:,event_idxs(i-1)+1:event_idxs(i));
            stm_t_eim_ei_eim = tmult(stm_t_eim_ei_0, stm_0_eim);
    
%             Qbart(:,:,event_idxs(i-1)+1:event_idxs(i)) = traj.Q_t(:,:,event_idxs(i-1)+1:event_idxs(i)) ...
%                 - tmult(  stm_t_eim_ei_eim  , tmult( Qbart(:,:,event_idxs(i-1)) , stm_t_eim_ei_eim, [0 1]  )    );

            Qbar_t(:,:,event_idxs(i-1)+1:event_idxs(i)) = traj.Q_t(:,:,event_idxs(i-1)+1:event_idxs(i)) ...
                - tmult(  stm_t_eim_ei_eim  , tmult( traj.Q_t(:,:,event_idxs(i-1)) , stm_t_eim_ei_eim, [0 1]  )    );


        end




    end

end

















%% old code for reference


% simparams.dynSys = dynSys;
% 
% 
% % Initial and final involved segment numbers
% segi = traj.t_s(idx_i);
% segf = traj.t_s(idx_f);
% 
% % Get the segment specific index for the first segment
% segi_idx_i = find_cell_t_i_idx_2(traj.t,traj.t_s,idx_i);
% 
% % If it is the last index in the specific segment, start at the beginning
% % of the next segment.
% fullInitialSeg = 0;
% if segi_idx_i == size(traj.stm_t_i{segi},3) && segf > segi
%     segi = segi + 1;
%     fullInitialSeg = 1;
% end
% 
% 
% 
% 
% % Initialize structure for Qbarf
% Qbarf = zeros(6,6);
% iter = 1;
% j = segf; % Start at the final segment and go backwards
% while iter
% 
% 
%     if j == segi && j == segf
%         % Get the effect from the beginning of segi/segf (segj) to idx_f
%         segj_idx_f = find_cell_t_i_idx_2(traj.t,traj.t_s,idx_f);
% 
%         Qbarf_f_segj0 = traj.Q_t_i{j}(:,:,segj_idx_f);
% 
%         % If don't start from the beginning of the segment, subtract the
%         % effect of Qbar on f from the beginning of the segment to
%         % segf_idx_i
% 
%         segj_idx_i = find_cell_t_i_idx_2(traj.t,traj.t_s,idx_i);
% 
%         if fullInitialSeg && segj_idx_i == size(traj.Q_t_i{segi-1},3)
%             segj_idx_i = 1;
%         end
% 
% 
%         if segj_idx_i > 1
%             % The Qbar effect at i from the beginning of seg j to i
%             Qbari_i_segj0 = traj.Q_t_i{segf}(:,:,segj_idx_i);
%             % STMs required to propagate the above effect to f
%             stm_i_j0 = traj.stm_t_i{j}(:,:,segj_idx_i);
%             stm_j0_i = invert_stm(stm_i_j0, simparams);
%             stm_f_j0 = traj.stm_t_i{j}(:,:,segj_idx_f);
%             stm_f_i = stm_f_j0 * stm_j0_i;
%             % The effect of Qbar at f from the beginning of segj to i
%             Qbarf_i_segj0 = stm_f_i * Qbari_i_segj0 * stm_f_i';
% 
%             % The total Qbar effect at f: subtracting the effect from j0 to
%             % i at f ... from the total effect from j0 to f at f
%             Qbarf = Qbarf_f_segj0 - Qbarf_i_segj0;
%         else
%             Qbarf = Qbarf_f_segj0;
%         end
% 
%         iter = 0;
%     elseif j == segf
%         % Include the effect from the beginning of segf to idx_f
%         segf_idx_f = find_cell_t_i_idx_2(traj.t,traj.t_s,idx_f);
% 
%         Qbarf_f_segf0 = traj.Q_t_i{segf}(:,:,segf_idx_f);
% 
%         Qbarf = Qbarf + Qbarf_f_segf0;
% 
%         % Get the STM from the beginning of seg f to idx_f. This starts the
%         % STM reverse combination. Cidx stands for current index.
%         stm_f_Cidx = traj.stm_t_i{segf}(:,:,segf_idx_f);
% 
%     elseif j > segi
%         % Include the entire effect of the full segment between segf & segi
%         Qbarj = traj.Q_i(:,:,j);
%         
%         Qbarf_jf_j0 = stm_f_Cidx * Qbarj * stm_f_Cidx';
% 
%         Qbarf = Qbarf + Qbarf_jf_j0;
% 
%         % Continue assembling the reverse total STM
% %         stm_jf_j0 = traj.stm_i(:,:,j);
%         stm_jf_j0 = traj.stm_t_i{j}(:,:,end);
% 
%         stm_f_Cidx = stm_f_Cidx * stm_jf_j0;
% 
% 
%     elseif j == segi
%         % Check if the whole segment or a partial segment contribution is
%         % being included
% 
%         if fullInitialSeg
%             % Include the entire effect of the full segment between segf &
%             % segi (same as above full segment)
%             Qbarj = traj.Q_i(:,:,j);
%             
%             Qbarf_jf_j0 = stm_f_Cidx * Qbarj * stm_f_Cidx';
%     
%             Qbarf = Qbarf + Qbarf_jf_j0;
%             
%         else
%             % Only include a fraction of the initial segment
% %             QbarjF_jf_i = QbarjF_jf_j0 - stm_jf_i * QbarjF_i_j0 * stm_jf_i'
%             % Get the segment index corresponding to idx_i
%             segj_idx_i = find_cell_t_i_idx_2(traj.t,traj.t_s,idx_i);
%             % Get the STM from idx_i to the end of the segment
%             stm_i_j0 = traj.stm_t_i{j}(:,:,segj_idx_i);
%             stm_j0_i = invert_stm(stm_i_j0, simparams);
% %             stm_jf_j0 = traj.stm_i(:,:,j);
%             stm_jf_j0 = traj.stm_t_i{j}(:,:,end);
%             stm_jf_i = stm_jf_j0 * stm_j0_i;
% 
%             % Get Qbar at the end of segment j (jf) for the entire segment
%             Qbarjf_jf_j0 = traj.Q_i(:,:,j);
% 
%             % Get the Qbar at i from the effect j0 to i
%             Qbari_i_j0 = traj.Q_t_i{j}(:,:,segj_idx_i);
%             % Propagate that to the end of segment j (jf)
%             Qbarjf_i_j0 = stm_jf_i * Qbari_i_j0 * stm_jf_i';
% 
%             % Subtract the impact from j0 to i at jf from the entire
%             % segment's impact at jf
%             Qbarjf_jf_i = Qbarjf_jf_j0 - Qbarjf_i_j0;
% 
%             % Propagate that to idx_f
%             Qbarf_jf_i = stm_f_Cidx * Qbarjf_jf_i * stm_f_Cidx';
% 
%             % Add it to Qbarf
%             Qbarf = Qbarf + Qbarf_jf_i;
%             
%             ppp=1;
%         end
% 
% 
% 
%         % End loop if at the initial segment
%         iter = 0;
%     else
%         assert(0,'Something went wrong.');
%     end
% 
% 
% 
% 
%     
% 
%     % Decrement j
%     j = j - 1;
% end


end