function optimHistoryMovie(history, outputFilepath, simparams)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

figure;
x=history.x;
writerObj = VideoWriter(outputFilepath);
writerObj.FrameRate = 10;
open(writerObj);

for i = 1:size(x,3)
%     x_opt, x_t, t_s, simparams, tcm_idx
    [~, ~, x_t, stm_t, t, t_s] = createStateStmHistory(x(:,:,i), simparams);
    [~,tcm_idx] = opt_multiple_tcm(x(:,:,i), t, t_s, stm_t, simparams);
    clf;
    plotMultiSegTraj(x(:,:,i), x_t, t_s, simparams, tcm_idx);
    frame = getframe(gcf);
    writeVideo(writerObj, frame);
end

close(writerObj);



end

