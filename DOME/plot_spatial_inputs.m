close all;
Render.time_plot = 0:45:Simulation.Tmax;
Render.all_time = 0:Simulation.deltaT:Simulation.Tmax;
Render.cmap_inputs = linspace2([0,0,0], [0,0.3,.9], 100)'; 
Render.color_rods = [1 0 0];
Render.sim_c = [0.3 0.3 0.3];
Render.exp_c = [1 0 0];


%% Plot the position of each agent at specified time instants (Time_plot)
%%%both in silico and in vivo. Also, it plots the comparison between real and simulated histograms

for i=1:length(Render.time_plot)

    %%%%%%%%%%%%%%%%%%%%%%%%% IN SILICO EXPERIMENT %%%%%%%%%%%%%%%%%%%%%%%%% 

    %get the position of the agents
    cur_ind = max([Render.time_plot(i)/Simulation.deltaT,2]);
    [~,indices_inWindow] = getInWindow(squeeze(xVec(cur_ind,:,:)), Simulation.arena);
    x_cur = squeeze(xVec(cur_ind,indices_inWindow,:));
    x_prev = squeeze(xVec(max([(Render.time_plot(i)/Simulation.deltaT)-1,1]),indices_inWindow,:));

    figure("Position",[0 0 500 300])
    %Plot the light projected on the sample
    if isfield(Environment,'Inputs') && isfield(Environment.Inputs,'Points')
        plotEnvField(Environment.Inputs.Points, Environment.Inputs.Values, Simulation.arena, Render.cmap_inputs)
    end
    %Plot the trajectories
    if Simulation.drawTraj; plotTrajectory(xVec, false, [0,0.7,0.9], Simulation.drawTraj); end
    %Plot the rods in the environment
    plotSwarmInit(x_cur,Render.time_plot(i),inf,inf,Simulation.arena,Simulation.arena,false,false,false, 'rod', [50,20], x_prev,[1 0 0]);
    %Plot the boundary of the simulation
    if isfield(Environment,'boundary'); plotBoundary(Environment.boundary); end
    %Save the image created in the output folder
    if outputDir
        saveas(gcf, fullfile(output_path, sprintf("spat_%d",Render.time_plot(i))))
        saveas(gcf, fullfile(output_path, sprintf("spat_%d",Render.time_plot(i))),'png')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%% IN VIVO EXPERIMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % get distribution wrt light intensity
    mask = detectObjects(data_folder, background_sub, brightness_thresh,Render.time_plot(i)+30);


    figure("Position",[500 100 500 500])
    x_vec = linspace(window(1),window(2),size(mask,2));
    y_vec = linspace(window(3),window(4),size(mask,1));
    box on    
    hold on
    colormap(Render.cmap_inputs)
    %Plot the inputs
    imagesc(x_vec,y_vec,flip(u'))
    %Plot the agents detected
    I=imagesc(x_vec,y_vec,cat(3,Render.color_rods(1).*mask,Render.color_rods(2).*mask,Render.color_rods(3).*mask));
    set(I, 'AlphaData', mask);
    %Resize the plot
    axis('equal')
    axis(window)
    xticks([])
    yticks([])
    if outputDir
        saveas(gcf, fullfile(output_path, sprintf("spat_exp_%d",Render.time_plot(i))))
        saveas(gcf, fullfile(output_path, sprintf("spat_exp_%d",Render.time_plot(i))),'png')
    end


    %%%%%%%%%%%%%%%%%%%%%%%%% HISTOGTAM COMPARISON %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Get density of experimental data
    [density_by_input_exp, bins, norm_slope_exp(i), ~, ~, ~,~, ~] = agentsDensityByInput(Environment.Inputs.Points, Environment.Inputs.Values, mask, window);
    [density_by_input_sim, ~, norm_slope_sim(i), ~, ~, ~,~, ~] = agentsDensityByInput(Environment.Inputs.Points, Environment.Inputs.Values, x_cur, window);


    
    figure("Position",[500+500 100+130 500 250])
    %Compute TVD between the distributions
    tvd = 0.5 * norm(density_by_input_exp-density_by_input_sim,1); % Total Variation Distance
    hold on
    %Generate the histograms
    b_exp = bar((bins(1:end-1)+bins(2:end))/2,density_by_input_exp, 1, FaceColor = 'b', FaceAlpha = 0.5);
    b_sim = bar((bins(1:end-1)+bins(2:end))/2,density_by_input_sim, 1, FaceColor = 'k', FaceAlpha = 0.4);
    legend({'REAL','SIMULATED'},'FontSize',14)
    xlabel('Input intensity','FontSize',14)
    ylabel('Density of agents','FontSize',14)
    yticks([0:0.25:1]);
    text(mean(bins),max(density_by_input_exp)*1.10,['TVD=',num2str(tvd,'%.2f')],'HorizontalAlignment','center','FontSize',14)
    ylim([0,max(density_by_input_exp)*1.15])
    xlim([-0.1,1.1])
    xticks(round(bins,2))
    box
    if outputDir
        saveas(gcf, fullfile(output_path, sprintf("hist_diff_%d",Render.time_plot(i))))
        saveas(gcf, fullfile(output_path, sprintf("hist_diff_%d",Render.time_plot(i))),'png')
    end




end


%%%%%%%%%%%%%%%%%%%%%%%%% Photo-accumulation Index comp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(Render.all_time)
    %Comput PhI at every time instant

    mask = detectObjects(data_folder, background_sub, brightness_thresh,Render.all_time(i)+30);
    cur_ind = max([Render.all_time(i)/Simulation.deltaT,2]);
    [~,indices_inWindow] = getInWindow(squeeze(xVec(cur_ind,:,:)), Simulation.arena);
    x_cur = squeeze(xVec(cur_ind,indices_inWindow,:));
    %Get density of experimental data
    [~, ~, norm_slope_exp(i), ~, ~, ~,~, ~] = agentsDensityByInput(Environment.Inputs.Points, Environment.Inputs.Values, mask, window);
    [~, ~, norm_slope_sim(i), ~, ~, ~,~, ~] = agentsDensityByInput(Environment.Inputs.Points, Environment.Inputs.Values, x_cur, window);

end

%plot the in silico and in vivo PhI to compare the settling time
figure("Position",[500 100+150 500 300])
plot(Render.all_time,norm_slope_sim,'LineWidth',2,'Color',Render.sim_c);
hold on;
plot(Render.all_time,norm_slope_exp,'LineWidth',2,'Color',Render.exp_c);
legend({'REAL','SIMULATED'},'FontSize',14)
xlabel('Time [s]','FontSize',14)
ylabel('PhI','FontSize',14)
xlim([0, Simulation.Tmax]);
if outputDir
        saveas(gcf, fullfile(output_path, 'slope_comp'))
        saveas(gcf, fullfile(output_path, 'slope_comp'),'png')
end


