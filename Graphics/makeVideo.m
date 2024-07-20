function [] = makeVideo(xVec,deltaT,Environment,Render,output_file)

% parse inputs
time_instants = size(xVec,1);

% print progress
fprintf('Generating video: ')
printProgress(0,time_instants)

% make video file
video = VideoWriter(output_file,'MPEG-4');
video.FrameRate = Render.frameRate;
open(video);

% init figure
figure
if isfield(Environment,'Inputs') && isfield(Environment.Inputs,'Points')
    plotEnvField(Environment.Inputs.Points, Environment.Inputs.Values, Render.window)
end
if Render.drawTraj; plotTrajectory(xVec, false, [0,0.7,0.9], Render.drawTraj); end
plotSwarmInit(squeeze(xVec(1,:,:)), 0, inf, inf, Render.window, [Render.window(2)-Render.window(1), Render.window(4)-Render.window(3)]/2, false, false, false, Render.agentShape, Render.agentSize, Render.agentsColor, squeeze(2*xVec(1,:,:)-xVec(2,:,:)));
if isfield(Environment,'boundary'); plotBoundary(Environment.boundary); end

% save frame
currFrame = getframe(gcf);
writeVideo(video,currFrame);

% for each time instant
for i=2:time_instants
    printProgress(i,time_instants)
    
    %clear figure
    cla
    t = (i-1)*deltaT;
    [x,indices_inWindow] = getInWindow(squeeze(xVec(i,:,:)), Render.window);
    x_pre=squeeze(xVec(i-1,indices_inWindow,:));
    
    % plot current positions and inputs
    if isfield(Environment,'Inputs') && isfield(Environment.Inputs,'Points')
        plotEnvField(Environment.Inputs.Points, Environment.Inputs.Values, Render.window)
    end
    if Render.drawTraj; plotTrajectory(xVec, false, [0,0.7,0.9], Render.drawTraj); end
    if isfield(Environment,'boundary'); plotBoundary(Environment.boundary); end
    plotSwarm(x, t, inf, inf, false, ones(size(x,1), 1), false, Render.agentShape, Render.agentSize, Render.agentsColor, x_pre);
    
    % save frame
    currFrame = getframe(gcf);
    writeVideo(video,currFrame);
end

% terminate video
close(video);

end

