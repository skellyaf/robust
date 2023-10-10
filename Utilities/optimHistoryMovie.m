function optimHistoryMovie(history, outputFilepath, simparams)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

figure;
x=history.x;
writerObj = VideoWriter(outputFilepath);
writerObj.FrameRate = 10;
open(writerObj);

nFrames = 50;

if size(x,3) < nFrames
    nFrames = size(x,3);

end
    

for i = round(linspace(1,size(x,3),nFrames))
% for i = 1:size(x,3)
%     x_opt, x_t, t_s, simparams, tcm_idx
%     [~, ~, x_t, stm_t, t, t_s] = createStateStmHistory(x(:,:,i), simparams);
    traj = createStateStmSttQdQHistory(x(:,:,i), simparams);

%     [~,tcm_idx] = opt_multiple_tcm(x(:,:,i), t, t_s, stm_t, simparams);
    [deltaV, deltaVs_nom] = calcDeltaV(x(:,:,i), traj.x_i_f, traj.stm_i, simparams);
    [~, tcm_idx] = opt_multiple_tcm_wQ(x(:,:,i), traj, deltaVs_nom, simparams);    
    clf;
    plotMultiSegTraj(x(:,:,i), traj.x_t, traj.t_s, simparams, tcm_idx);
    frame = getframe(gcf);
    writeVideo(writerObj, frame);
end

close(writerObj);



end

