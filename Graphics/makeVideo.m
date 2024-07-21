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
fig.Position = [100 100 2000 1500];

[x,indices_inWindow] = getInWindow(squeeze(xVec(1,:,:)), Render.window);
x_pre = squeeze(2*xVec(1,indices_inWindow,:)-xVec(2,indices_inWindow,:));

% plot swarm
plotSwarmInit(x, 0, inf, inf, Render.window, [Render.window(2)-Render.window(1), Render.window(4)-Render.window(3)]/2, false, false, false, Render.agentShape, Render.agentSize, Render.agentsColor, x_pre);
makeFigureSwarm(x, x_pre, 0,Environment, Render)

% plot light distribution
% makeFigureDensity(x,0,Environment,Render)

% plot experimental image
% makeFigureExperiment(0, Render)

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
    
    plot current positions and inputs
    makeFigureSwarm(x, x_pre, t,Environment, Render)
    
    % plot light distribution
    % makeFigureDensity(x,t,Environment,Render)
    
    % plot experimental image
    % makeFigureExperiment(t, Render)
    
    % save frame
    currFrame = getframe(gcf);
    writeVideo(video,currFrame);
    
end

% terminate video
close(video);

end

%% Functions to print a single frame
% print frame from simulation
function [] = makeFigureSwarm(x, x_pre, t,Environment, Render)
if isfield(Environment,'Inputs') && isfield(Environment.Inputs,'Points')
    plotEnvField(Environment.Inputs.Points, Environment.Inputs.Values, Render.window, Render.cmap_inputs)
end
if Render.drawTraj; plotTrajectory(xVec, false, [0,0.7,0.9], Render.drawTraj); end
if isfield(Environment,'boundary'); plotBoundary(Environment.boundary); end
plotSwarm(x, t, inf, inf, false, ones(size(x,1), 1), false, Render.agentShape, Render.agentSize, Render.agentsColor, x_pre);
end

% print light distribution from simulation
function [] = makeFigureDensity(x, t, Environment, Render)
[density_by_input_sim, bins] = agentsDensityByInput(Environment.Inputs.Points, Environment.Inputs.Values, x, Render.window, 5);

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

% print frame from experiment
function [] = makeFigureExperiment(t, Render)
experiments_folder = "/Volumes/DOMEPEN/Experiments";
experiments_names = ["2023_07_10_E_30","2023_07_10_E_34","2023_07_10_E_35"];
brightness_thresh = 0.3;    % brightness threshold for object detection (spatial experiments)
background_sub = true;      % use background subtraction for object detection (spatial experiments)
pattern_blurring = 15;      % blurring of the spatial input pattern (spatial experiments)
exp_setup_time = 30;     % initial time window to be discarded from the experiment [s] (set to 30 for BCL, 0 for other exp)

mask = cell(size(experiments_names));
u = cell(size(experiments_names));

% combo exp data
for j=1:length(experiments_names)   % for each replicate
    experiments_names(j) = strrep(experiments_names(j),'_E_','_Euglena_');
    data_folder =  fullfile(experiments_folder,experiments_names(j));
    
    time_to_plot_exp = t+exp_setup_time;
    mask{j} = detectObjects(data_folder, background_sub, brightness_thresh, time_to_plot_exp);
    u{j} = loadInputPattern(data_folder, pattern_blurring, time_to_plot_exp);
    assert( mean(abs(u{1} - imresize(u{j},size(u{1}))),'all')<0.1, 'Replicates have different inputs' )
    
    % get distribution wrt light intensity
    %         [~, ~, agents_by_input(j,:), pixels_by_input(j,:)] = agentsDensityByInput(inputs.Points, inputs.Values, mask{j}, Render.window, n_bins);
    
    % combination of masks over the replicates
    if j==1
        combo_mask = mask{1};
    else
        combo_mask = (combo_mask+mask{j})>=1;
    end
end

%     % weighted average light distribution over the replicates
%     density_by_input_exp_mean = sum(squeeze(agents_by_input)./squeeze(pixels_by_input),1);
%     density_by_input_exp_mean = density_by_input_exp_mean/sum(density_by_input_exp_mean);
%     [~,norm_slope_exp(i),~] = linearDependence((bins(1:end-1)+bins(2:end))'/2, density_by_input_exp_mean');

x_vec = linspace(Render.window(1),Render.window(2),size(combo_mask,2));
y_vec = linspace(Render.window(3),Render.window(4),size(combo_mask,1));
%Plot the inputs
imagesc(x_vec,y_vec,flip(u{1}'))
%Plot the agents detected
I=imagesc(x_vec,y_vec,cat(3,Render.agentsColor(1).*combo_mask,Render.agentsColor(2).*combo_mask,Render.agentsColor(3).*combo_mask));
set(I, 'AlphaData', combo_mask);
title(sprintf('t=%.2fs',time_to_plot_exp))
end
