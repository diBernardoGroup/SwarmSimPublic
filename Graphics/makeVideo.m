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
fig=figure;
fig.Color = 'w';
[x,indices_inWindow] = getInWindow(squeeze(xVec(1,:,:)), Render.window);
x_pre = squeeze(2*xVec(1,indices_inWindow,:)-xVec(2,indices_inWindow,:));
    
% plot swarm
plotSwarmInit(x, 0, inf, inf, Render.window, [Render.window(2)-Render.window(1), Render.window(4)-Render.window(3)]/2, false, false, false, Render.agentShape, Render.agentSize, Render.agentsColor, x_pre);
makeFigureSwarm(x, x_pre, 0,Environment, Render)

% plot light distribution
% makeFigureDensity(x,0,Environment,Render)

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
     makeFigureSwarm(x, x_pre, t,Environment, Render)
    
    % plot light distribution
    % makeFigureDensity(x,t,Environment,Render)
    
    % save frame
    currFrame = getframe(gcf);
    writeVideo(video,currFrame);
    
end

% terminate video
close(video);

end

function [] = makeFigureSwarm(x, x_pre, t,Environment, Render)
    if isfield(Environment,'Inputs') && isfield(Environment.Inputs,'Points')
        plotEnvField(Environment.Inputs.Points, Environment.Inputs.Values, Render.window, Render.cmap_inputs)
    end
    if Render.drawTraj; plotTrajectory(xVec, false, [0,0.7,0.9], Render.drawTraj); end
    if isfield(Environment,'boundary'); plotBoundary(Environment.boundary); end
    plotSwarm(x, t, inf, inf, false, ones(size(x,1), 1), false, Render.agentShape, Render.agentSize, Render.agentsColor, x_pre);
end

function [] = makeFigureDensity(x, t, Environment, Render)
    % simulation light distribution
    [density_by_input_sim, bins, norm_slope_sim, c_coeff_sim, coefficents, ~,~, u_values_sim] = agentsDensityByInput(Environment.Inputs.Points, Environment.Inputs.Values, x, Render.window, 5);

    title(sprintf('t=%.2fs',t))
    bar((bins(1:end-1)+bins(2:end))/2,density_by_input_sim, 1, 'FaceColor', [1 1 1]*0.5)
    hold on
    plot(bins,coefficents(1)+coefficents(2)*bins,LineWidth=2);
    xlabel('Input intensity')
    ylabel('Density of agents')
    yticks([0:0.2:1]);
    %text(max(bins),max(density_by_input_sim)*1.1,['\rho=',num2str(c_coeff_sim,'%.2f')],'HorizontalAlignment','right','FontSize',14)
    text(max(bins),0.35,['PhAI=',num2str(norm_slope_sim,'%.2f')],'HorizontalAlignment','right','FontSize',14)
    %ylim([0,max(density_by_input_sim)*1.15])
    ylim([0,0.4])
    xlim([-0.1,1.1])
    xticks(round(bins,2))
end


