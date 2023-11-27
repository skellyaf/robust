function optimHistoryMovie(history, outputFilepath, simparams)
%optimHistoryMovie Makes a movie (.avi) of the trajectory optimization
%convergence steps taken from the initial guess to the converged
%trajectory. 

% The history data structure includes the multiple segment trajectory
% parameter vector (x) for each iteration.

% outputFilePath is the path where you want the video output. Can just do
% outputFilepath = './' for current directory.

% Open a figure
figure;
% Intermediate variable
x=history.x;

% The matlab video output file creation
writerObj = VideoWriter(outputFilepath);

% I adjusted the framerate to pause a little longer - not required pending
% how many trajectory increments you plot
writerObj.FrameRate = 10;

% Open the video writer
open(writerObj);

% The entire history is too many frames in some cases, I plot every 50.
nFrames = 50;

% unless there are less than 50 frames
if size(x,3) < nFrames
    nFrames = size(x,3);
end
    
% Loop for each plot to add to the video
for i = round(linspace(1,size(x,3),nFrames))

    % Propagating the position vectors (and other stuff) for each x in the
    % history tensor (the x(:,:,i) part)
    traj = createStateStmSttQdQHistory(x(:,:,i), simparams);

    % Calculating nominal delta V and the optimal TCMs to plot / display
    [deltaV, deltaVs_nom] = calcDeltaV(x(:,:,i), traj.x_i_f, traj.stm_i, simparams);
    [~, tcm_idx] = opt_multiple_tcm_wQ(x(:,:,i), traj, deltaVs_nom, simparams);    

    % Clear the previous stuff that was on the figure
    clf;

    % Function that plots the trajectory
    plotMultiSegTraj(x(:,:,i), traj.x_t, traj.t_s, simparams, tcm_idx);

    % Getting the frame from the current plot
    frame = getframe(gcf);

    % Writing the frame to the video
    writeVideo(writerObj, frame);
end

% Closing the video writer object
close(writerObj);



end

