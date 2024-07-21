close all;
clear

defaultParamMicroorg;

sim_folder = '/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/media/model validation/Euglena/simulations/2024_06_25_BCLx36_1';

experiments_folder = "/Volumes/DOMEPEN/Experiments";
experiments_names = ["2023_07_10_E_30","2023_07_10_E_34","2023_07_10_E_35"];

sim_data = load(fullfile(sim_folder,'data.mat'));
Simulation = sim_data.Simulation;
Environment = sim_data.Environment;
inputs = sim_data.Environment.Inputs;
xVec = sim_data.xVec;
data_folder = sim_data.data_folder;

make_videos = true;
% Render.time_plot = 0:45:Simulation.Tmax;
Render.time_plot = [];
% Render.all_time  = 0:Simulation.deltaT:Simulation.Tmax;
Render.all_time  = 0:10:Simulation.Tmax;
Render.window = [-Simulation.arena(1),Simulation.arena(1),-Simulation.arena(2),Simulation.arena(2)]/2; % size of the simulation window

exp_setup_time = 30;     % initial time window to be discarded from the experiment [s] (set to 30 for BCL, 0 for other exp)
n_bins         = 5;      % number of bins for light distribution (set to 5 for BCL, 3 for other exp)

output_path = sim_folder;

%% Plot at selected time instants (Render.time_plot)
%%%both in silico and in vivo. Also, it plots the comparison between real and simulated histograms

for i=1:length(Render.time_plot)
    
    %%%%%%%%%%%%%%%%%%%%%%%%% IN SILICO EXPERIMENT %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %get the position of the agents
    frame_index = round(Render.time_plot(i) / Simulation.deltaT) + 1;
    [x_cur,indices_inWindow] = getInWindow(squeeze(xVec(frame_index,:,:)), Render.window);
    if frame_index == 1
        x_prev = squeeze(2*xVec(1,indices_inWindow,:)-xVec(2,indices_inWindow,:));
    else
        x_prev = squeeze(xVec(frame_index-1,indices_inWindow,:));
    end
    [density_by_input_sim, bins] = agentsDensityByInput(Environment.Inputs.Points, Environment.Inputs.Values, x_cur, Render.window, n_bins);
    
    figure("Position",[0 0 500 300])
    %Plot the light projected on the sample
    if isfield(Environment,'Inputs') && isfield(Environment.Inputs,'Points')
        plotEnvField(Environment.Inputs.Points, Environment.Inputs.Values, Render.window, Render.cmap_inputs)
    end
    %Plot the trajectories
    if Render.drawTraj; plotTrajectory(xVec, false, [0,0.7,0.9], Render.drawTraj); end
    %Plot the rods in the environment
    plotSwarmInit(x_cur,Render.time_plot(i),inf,inf,Render.window,Simulation.arena,false,false,false, Render.agentShape, Render.agentSize, Render.agentsColor, x_prev);
    %Plot the boundary of the simulation
    if isfield(Environment,'boundary'); plotBoundary(Environment.boundary); end
    %Save the image created in the output folder
    if output_path
        saveas(gcf, fullfile(output_path, sprintf("positions_sim_%d",Render.time_plot(i))))
        saveas(gcf, fullfile(output_path, sprintf("positions_sim_%d",Render.time_plot(i))),'png')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%% IN VIVO EXPERIMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % get distribution wrt light intensity
    %     mask = detectObjects(data_folder, background_sub, brightness_thresh,Render.time_plot(i)+exp_setup_time);
    
    for j=1:length(experiments_names)   % for each replicate
        time_to_plot_exp = Render.time_plot(i)+exp_setup_time;
        experiments_names(j) = strrep(experiments_names(j),'_E_','_Euglena_');
        data_folder =  fullfile(experiments_folder,experiments_names(j));
        mask{j} = detectObjects(data_folder, background_sub, brightness_thresh, time_to_plot_exp);
        u{j} = loadInputPattern(data_folder, pattern_blurring, time_to_plot_exp);
        assert( mean(abs(u{1} - imresize(u{j},size(u{1}))),'all')<0.1, 'Replicates have different inputs' )
        
        % get distribution wrt light intensity
        [~, ~, agents_by_input(j,:), pixels_by_input(j,:)] = agentsDensityByInput(inputs.Points, inputs.Values, mask{j}, Render.window, n_bins);
        
        % combination of masks over the replicates
        if j==1
            combo_mask = mask{1};
        else
            combo_mask = (combo_mask+mask{j})>=1;
        end
    end
    
    % weighted average light distribution over the replicates
    density_by_input_exp_mean = sum(squeeze(agents_by_input)./squeeze(pixels_by_input),1);
    density_by_input_exp_mean = density_by_input_exp_mean/sum(density_by_input_exp_mean);
    
    
    figure("Position",[500 100 500 500])
    x_vec = linspace(Render.window(1),Render.window(2),size(combo_mask,2));
    y_vec = linspace(Render.window(3),Render.window(4),size(combo_mask,1));
    box on
    hold on
    colormap(Render.cmap_inputs)
    %Plot the inputs
    imagesc(x_vec,y_vec,flip(u{1}'))
    %Plot the agents detected
    I=imagesc(x_vec,y_vec,cat(3,Render.agentsColor(1).*combo_mask,Render.agentsColor(2).*combo_mask,Render.agentsColor(3).*combo_mask));
    set(I, 'AlphaData', combo_mask);
    %Resize the plot
    axis('equal')
    axis(Render.window)
    xticks([])
    yticks([])
    if output_path
        saveas(gcf, fullfile(output_path, sprintf("positions_exp_%d",Render.time_plot(i))))
        saveas(gcf, fullfile(output_path, sprintf("positions_exp_%d",Render.time_plot(i))),'png')
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%% HISTOGTAM COMPARISON %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure("Position",[500+500 100+130 500 250])
    %Compute TVD between the distributions
    tvd = 0.5 * norm(density_by_input_exp_mean-density_by_input_sim,1); % Total Variation Distance
    hold on
    %Generate the histograms
    b_exp = bar((bins(1:end-1)+bins(2:end))/2,density_by_input_exp_mean, 1, FaceColor = 'b', FaceAlpha = 0.5);
    b_sim = bar((bins(1:end-1)+bins(2:end))/2,density_by_input_sim, 1, FaceColor = 'k', FaceAlpha = 0.4);
    legend({'REAL','SIMULATED'},'FontSize',14)
    xlabel('Input intensity','FontSize',14)
    ylabel('Density of agents','FontSize',14)
    yticks([0:0.25:1]);
    text(mean(bins),max(density_by_input_exp_mean)*1.10,['TVD=',num2str(tvd,'%.2f')],'HorizontalAlignment','center','FontSize',14)
    ylim([0,max(density_by_input_exp_mean)*1.15])
    xlim([-0.1,1.1])
    xticks(round(bins,2))
    box
    if output_path
        saveas(gcf, fullfile(output_path, sprintf("difference_light_distribution_%d",Render.time_plot(i))))
        saveas(gcf, fullfile(output_path, sprintf("difference_light_distribution_%d",Render.time_plot(i))),'png')
    end
    
    
    
    
end

%% Plot at each time instant (Render.all_time)

% plot PhAI in time
fig_phai=figure("Position",[500 100 500 300]);
fig_phai.Color = 'w';
hold on;
xlabel('Time [s]','FontSize',14)
ylabel('PhAI','FontSize',14)
xlim([0, Simulation.Tmax]);
ylim([-0.6,0.6])
yticks([-0.6:0.2:0.6])
box on
grid on
% make video file
if make_videos
    video_phai = VideoWriter(fullfile(output_path, 'video_phai'),'MPEG-4');
    video_phai.FrameRate = Render.frameRate;
    open(video_phai);
end

% plot TVD in time
fig_tvd=figure("Position",[500 100+300 500 300]);
fig_tvd.Color = 'w';
hold on;
xlabel('Time [s]','FontSize',14)
ylabel('TVD','FontSize',14)
xlim([0, Simulation.Tmax]);
ylim([0,0.15])
box on
grid on
% make video file
if make_videos
    video_tvd = VideoWriter(fullfile(output_path, 'video_tvd'),'MPEG-4');
    video_tvd.FrameRate = Render.frameRate;
    open(video_tvd);
end

fig_light_dist=figure("Position",[500 100+300*2 500 300]);
fig_light_dist.Color = 'w';
hold on;
xlabel('Input intensity','FontSize',14)
ylabel('Density of agents','FontSize',14)
ylim([0,0.4])
xlim([-0.1,1.1])
xticks(0:1/n_bins:1)
box
% make video file
if make_videos
    video_light_dist = VideoWriter(fullfile(output_path, 'video_diff_light_dist'),'MPEG-4');
    video_light_dist.FrameRate = Render.frameRate;
    open(video_light_dist);
end

% allocate vars
norm_slope_exp = nan(size(Render.all_time));
norm_slope_sim = nan(size(Render.all_time));
tvd            = nan(size(Render.all_time));

for i=1:length(Render.all_time)
    % load simulation data
    frame_index = round(Render.all_time(i) / Simulation.deltaT) + 1;
    [x_cur,indices_inWindow] = getInWindow(squeeze(xVec(frame_index,:,:)), Render.window);
    
    [density_by_input_sim,bins] = agentsDensityByInput(Environment.Inputs.Points, Environment.Inputs.Values, x_cur, Render.window, n_bins);
    [~,norm_slope_sim(i),~] = linearDependence((bins(1:end-1)+bins(2:end))'/2, density_by_input_sim');
    
    % load and combo exp data
    for j=1:length(experiments_names)   % for each replicate
        time_to_plot_exp = Render.all_time(i)+exp_setup_time;
        experiments_names(j) = strrep(experiments_names(j),'_E_','_Euglena_');
        data_folder =  fullfile(experiments_folder,experiments_names(j));
        mask{j} = detectObjects(data_folder, background_sub, brightness_thresh, time_to_plot_exp);
        u{j} = loadInputPattern(data_folder, pattern_blurring, time_to_plot_exp);
        assert( mean(abs(u{1} - imresize(u{j},size(u{1}))),'all')<0.1, 'Replicates have different inputs' )
        
        % get distribution wrt light intensity
        [~, ~, agents_by_input(j,:), pixels_by_input(j,:)] = agentsDensityByInput(inputs.Points, inputs.Values, mask{j}, Render.window, n_bins);
        
        % combination of masks over the replicates
        if j==1
            combo_mask = mask{1};
        else
            combo_mask = (combo_mask+mask{j})>=1;
        end
    end
    
    % weighted average light distribution over the replicates
    density_by_input_exp_mean = sum(squeeze(agents_by_input)./squeeze(pixels_by_input),1);
    density_by_input_exp_mean = density_by_input_exp_mean/sum(density_by_input_exp_mean);
    [~,norm_slope_exp(i),~] = linearDependence((bins(1:end-1)+bins(2:end))'/2, density_by_input_exp_mean');
    
    %Get density of experimental data
    tvd(i) = 0.5 * norm(density_by_input_exp_mean-density_by_input_sim,1); % Total Variation Distance
    
    
    if make_videos || i==length(Render.all_time)
        
        % plot the in silico and in vivo PhAI to compare the settling time
        figure(fig_phai)
        cla
        plot(Render.all_time,norm_slope_exp,'LineWidth',2,'Color',Render.exp_c);
        plot(Render.all_time,norm_slope_sim,'LineWidth',2,'Color',Render.sim_c);
        % legend({'REAL','SIMULATED'},'FontSize',14)
        legend({'In vivo','In silico'},'FontSize',14,'FontAngle','italic')
        % save frame
        if make_videos
            currFrame = getframe(fig_phai);
            writeVideo(video_phai,currFrame);
        end
        
        % plot TVD in time
        figure(fig_tvd)
        cla
        plot(Render.all_time,tvd,'LineWidth',2,'Color', [1 1 1]*0.5);
        % save frame
        if make_videos
            currFrame = getframe(fig_tvd);
            writeVideo(video_tvd,currFrame);
        end
        
        % plot histograms distributions with respect to light intensity
        figure(fig_light_dist)
        cla
        b_exp = bar((bins(1:end-1)+bins(2:end))/2,density_by_input_exp_mean, 1, FaceColor = 'b', FaceAlpha = 0.5);
        b_sim = bar((bins(1:end-1)+bins(2:end))/2,density_by_input_sim, 1, FaceColor = 'k', FaceAlpha = 0.4);
        %legend({'REAL','SIMULATED'},'FontSize',14)
        legend({'In vivo','In silico'},'FontSize',14,'FontAngle','italic')
        text(mean(bins),0.35,['TVD=',num2str(tvd(i),'%.2f')],'HorizontalAlignment','center','FontSize',14)
        % save frame
        if make_videos
            currFrame = getframe(fig_light_dist);
            writeVideo(video_light_dist,currFrame);
        end
        
    end
end

% terminate videos
if make_videos
    close(video_phai);
    close(video_tvd);
    close(video_light_dist);
end

% save figures at the last time step
if output_path
    fig_phai.Units = fig_phai.PaperUnits; fig_phai.PaperSize = fig_phai.Position(3:4); % set correct pdf size
    saveas(fig_phai, fullfile(output_path, "photo_acc_index"))
    saveas(fig_phai, fullfile(output_path, "photo_acc_index"),'pdf')
    
    fig_tvd.Units = fig_tvd.PaperUnits; fig_tvd.PaperSize = fig_tvd.Position(3:4); % set correct pdf size
    saveas(fig_tvd, fullfile(output_path, 'TVD'))
    saveas(fig_tvd, fullfile(output_path, 'TVD'),'pdf')
    
    fig_light_dist.Units = fig_light_dist.PaperUnits; fig_light_dist.PaperSize = fig_light_dist.Position(3:4); % set correct pdf size
    saveas(fig_light_dist, fullfile(output_path, 'difference_light_distribution'))
    saveas(fig_light_dist, fullfile(output_path, 'difference_light_distribution'),'pdf')
end
