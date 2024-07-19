function [] = makeVideo(xVec,Environment,Render,Simulation,output_file)

time_instants = size(xVec,1);

% make video file
video = VideoWriter(output_file,'MPEG-4');
video.FrameRate = Render.frameRate;
open(video);

% init figure
figure
if isfield(Environment,'Inputs') && isfield(Environment.Inputs,'Points')
    x_vec = linspace(Render.window(1),Render.window(2),1000);
    y_vec = linspace(Render.window(3),Render.window(4),1000);
    [x_mesh, y_mesh] = meshgrid(x_vec, y_vec);
    F = griddedInterpolant(Environment.Inputs.Points,Environment.Inputs.Values, 'linear', 'nearest');
    imagesc(x_vec,y_vec,F(x_mesh',y_mesh')')
    colormap(Render.cmap_inputs)
end
plotSwarmInit(squeeze(xVec(1,:,:)), 0, inf, inf, Render.window, [Render.window(2)-Render.window(1), Render.window(4)-Render.window(3)]/2, false, false, false, Render.agentShape, Render.agentSize, Render.agentsColor, squeeze(2*xVec(1,:,:)-xVec(2,:,:)));
% save frame
currFrame = getframe(gcf);
writeVideo(video,currFrame);

fprintf('Generating video: ')
printProgress(0,time_instants)

% for each time instant
for i=2:time_instants
    printProgress(i,time_instants)
    
    %clear figure
    cla
    t = (i-1)*Simulation.deltaT;
    x=squeeze(xVec(i,:,:));
    x_pre=squeeze(xVec(i-1,:,:));
    
    % plot current positions and inputs
    if isfield(Environment,'Inputs')
        if isfield(Environment.Inputs,'Points')
            imagesc(x_vec,y_vec,F(x_mesh',y_mesh')')
        elseif isfield(Environment.Inputs,'Times')
            set(gca,'Color',interp1(linspace(0,1,100),Render.cmap_inputs,interp1(Environment.Inputs.Times, Environment.Inputs.Values, t, Environment.Inputs.InterpMethod)))
        end
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

