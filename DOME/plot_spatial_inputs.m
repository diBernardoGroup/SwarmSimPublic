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
% u_sim = sim_data.u;

Render.time_plot = 0:180:Simulation.Tmax;
Render.all_time  = 0:1:Simulation.Tmax;
Render.window = [-Simulation.arena(1),Simulation.arena(1),-Simulation.arena(2),Simulation.arena(2)]/2; % size of the simulation window

exp_setup_time = 30;     % initial time window to be discarded from the experiment [s] (set to 30 for BCL, 0 for other exp)
n_bins         = 5;     % number of bins for light distribution


% Render.window = [-Simulation.arena(1),Simulation.arena(1),-Simulation.arena(2),Simulation.arena(2)]/2; % size of the simulation window
% Render.drawON=false;                % draw swarm during simulation (if N is large slows down the simulation)
% Render.drawTraj=0;                  % draw trajectories of the agents (if N is large slows down the simulation)
% Render.recordVideo=false;           % record video of the simulation (if true drawON must be true)
% Render.frameRate = 1/Simulation.deltaT;
% Render.agentShape = "rod";          % shape to plot the agents "rod" or any defualt marker key ('.','+','diamond',...)
% Render.agentSize = 30;
% Render.all_time = 0:Simulation.deltaT:Simulation.Tmax;
% Render.shaded = true;
%
% % Light palette - red inputs
% Render.agentsColor = [0 0 1]; % blue
% Render.sim_c =      [0 0 0]; % black
% Render.exp_c =      [0 0 1]; % blue
% Render.cmap_inputs = linspace2([1,1,1], [1,0.5,0.5], 100)';  % light red

output_path = sim_folder;

%% Plot the position of each agent at specified time instants (Time_plot)
%%%both in silico and in vivo. Also, it plots the comparison between real and simulated histograms

for i=1:length(Render.time_plot)
    
    %%%%%%%%%%%%%%%%%%%%%%%%% IN SILICO EXPERIMENT %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %get the position of the agents
    cur_ind = max([Render.time_plot(i)/Simulation.deltaT,2]);
    [~,indices_inWindow] = getInWindow(squeeze(xVec(cur_ind,:,:)), Simulation.arena);
    x_cur = squeeze(xVec(cur_ind,indices_inWindow,:));
    x_prev = squeeze(xVec(max([(Render.time_plot(i)/Simulation.deltaT)-1,1]),indices_inWindow,:));
    
    [density_by_input_sim, bins, norm_slope_sim(i), ~, ~, ~,~, ~] = agentsDensityByInput(Environment.Inputs.Points, Environment.Inputs.Values, x_cur, Render.window, n_bins);
    
    %     figure("Position",[0 0 500 300])
    %     %Plot the light projected on the sample
    %     if isfield(Environment,'Inputs') && isfield(Environment.Inputs,'Points')
    %         plotEnvField(Environment.Inputs.Points, Environment.Inputs.Values, Simulation.arena, Render.cmap_inputs)
    %     end
    %     %Plot the trajectories
    %     if Render.drawTraj; plotTrajectory(xVec, false, [0,0.7,0.9], Render.drawTraj); end
    %     %Plot the rods in the environment
    %     plotSwarmInit(x_cur,Render.time_plot(i),inf,inf,Render.window,Simulation.arena,false,false,false, Render.agentShape, Render.agentSize, Render.agentsColor, x_prev);
    %     %Plot the boundary of the simulation
    %     if isfield(Environment,'boundary'); plotBoundary(Environment.boundary); end
    %     %Save the image created in the output folder
    %     if output_path
    %         saveas(gcf, fullfile(output_path, sprintf("spat_%d",Render.time_plot(i))))
    %         saveas(gcf, fullfile(output_path, sprintf("spat_%d",Render.time_plot(i))),'png')
    %     end
    
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
        [~, ~, ~, ~, ~, agents_by_input(j,:), pixels_by_input(j,:)] = agentsDensityByInput(inputs.Points, inputs.Values, mask{j}, Render.window, n_bins);
        
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
        saveas(gcf, fullfile(output_path, sprintf("spat_exp_%d",Render.time_plot(i))))
        saveas(gcf, fullfile(output_path, sprintf("spat_exp_%d",Render.time_plot(i))),'png')
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%% HISTOGTAM COMPARISON %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Get density of experimental data
    %[density_by_input_exp, bins, norm_slope_exp(i), ~, ~, ~,~, ~] = agentsDensityByInput(Environment.Inputs.Points, Environment.Inputs.Values, mask, Render.window);
    
    
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
        saveas(gcf, fullfile(output_path, sprintf("hist_diff_%d",Render.time_plot(i))))
        saveas(gcf, fullfile(output_path, sprintf("hist_diff_%d",Render.time_plot(i))),'png')
    end
    
    
    
    
end


%%%%%%%%%%%%%%%%%%%%%%%%% Photo-accumulation Index comp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(Render.all_time)
    %Comput PhI at every time instant
    
    % simulation data
    cur_ind = max([Render.all_time(i)/Simulation.deltaT,2]);
    [~,indices_inWindow] = getInWindow(squeeze(xVec(cur_ind,:,:)), Simulation.arena);
    x_cur = squeeze(xVec(cur_ind,indices_inWindow,:));
    [density_by_input_sim, ~, norm_slope_sim(i), ~, ~, ~,~, ~] = agentsDensityByInput(Environment.Inputs.Points, Environment.Inputs.Values, x_cur, Render.window, n_bins);
    
    % combo exp data
    for j=1:length(experiments_names)   % for each replicate
        time_to_plot_exp = Render.all_time(i)+exp_setup_time;
        experiments_names(j) = strrep(experiments_names(j),'_E_','_Euglena_');
        data_folder =  fullfile(experiments_folder,experiments_names(j));
        mask{j} = detectObjects(data_folder, background_sub, brightness_thresh, time_to_plot_exp);
        u{j} = loadInputPattern(data_folder, pattern_blurring, time_to_plot_exp);
        assert( mean(abs(u{1} - imresize(u{j},size(u{1}))),'all')<0.1, 'Replicates have different inputs' )
        
        % get distribution wrt light intensity
        [~, ~, ~, ~, ~, agents_by_input(j,:), pixels_by_input(j,:)] = agentsDensityByInput(inputs.Points, inputs.Values, mask{j}, Render.window, n_bins);
        
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
    
    
    %Get density of experimental data
    tvd(i) = 0.5 * norm(density_by_input_exp_mean-density_by_input_sim,1); % Total Variation Distance
end

% %plot the in silico and in vivo PhAI to compare the settling time
% figure("Position",[500 100+150 500 300])
% hold on;
% plot(Render.all_time,norm_slope_exp,'LineWidth',2,'Color',Render.exp_c);
% plot(Render.all_time,norm_slope_sim,'LineWidth',2,'Color',Render.sim_c);
% legend({'REAL','SIMULATED'},'FontSize',14)
% xlabel('Time [s]','FontSize',14)
% ylabel('PhAI','FontSize',14)
% xlim([0, Simulation.Tmax]);
% box on
% grid on

figure("Position",[500 100+450 500 300])
plot(Render.all_time,tvd,'LineWidth',2,'Color', [1 1 1]*0.5);
xlabel('Time [s]','FontSize',14)
ylabel('TVD','FontSize',14)
xlim([0, Simulation.Tmax]);
ylim([0, 0.12])
box on
grid on
if output_path
    fig=gcf; fig.Units = fig.PaperUnits; fig.PaperSize = fig.Position(3:4); % set correct pdf size
    saveas(gcf, fullfile(output_path, 'slope_comp'))
    saveas(gcf, fullfile(output_path, 'slope_comp'),'png')
    saveas(gcf, fullfile(output_path, 'slope_comp'),'pdf')
end


