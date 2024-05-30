clear
close all

defaultParamMicroorg

simulations_folder = '/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations';
experiments_folder = "/Volumes/DOMEPEN/Experiments";


% half_half
tags = ["half_half"];
sim_names = ["2024_05_30_half_half_1"];
experiments_names = ["2023_06_12_Euglena_2","2023_06_14_Euglena_6","2023_06_15_Euglena_12"];%,"2023_06_26_Euglena_29","2023_06_26_Euglena_30","2023_06_23_Euglena_1","2023_06_23_Euglena_2","2023_06_26_Euglena_2","2023_06_26_Euglena_1"];

% spatial
tags = ["half_half","grad_centr_light","grad_centr_dark"];
sim_names = ["2024_05_30_half_half_1";"2024_05_30_grad_centr_light_2";"2024_05_30_grad_centr_dark_1"];
experiments_names = ["2023_06_12_E_2","2023_06_14_E_6","2023_06_15_E_12";
    "2023_06_12_E_3","2023_06_12_E_4","2023_06_14_E_7";
    "2023_06_14_E_10","2023_06_15_E_15","2023_06_23_E_7"];
output_folder = '/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/comparison/Euglena spatial';

deltaT = 0.5;
timeInstants = [0:deltaT:180];


%% LOAD DATA
combo_mask = cell(1,size(experiments_names,1));
assert((length(tags)==length(sim_names)) && (length(sim_names)==size(experiments_names,1)))
for i = 1:size(experiments_names,1)  % for each experiment
    
    % load simulation data
    sim_folder = fullfile(simulations_folder,sim_names(i));
    sim_data = load(fullfile(sim_folder,'data.mat'));
    xFinal_inWindow{i} = sim_data.xFinal_inWindow;
    arena = sim_data.Simulation.arena;
    window = [-arena(1),arena(1),-arena(2),arena(2)]/2;
    inputs{i} = sim_data.Environment.Inputs;
    [density_by_input_sim(i,:), bins, norm_slope_sim(i), c_coeff_sim(i)] = agentsDensityByInput(inputs{i}.Points, inputs{i}.Values, xFinal_inWindow{i}, window);
    
    % load experiment data
    for j=1:size(experiments_names,2)   % for each replicate
        experiments_names(i,j) = strrep(experiments_names(i,j),'_E_','_Euglena_');
        data_folder =  fullfile(experiments_folder,experiments_names(i,j));
        mask{i,j} = detectObjects(data_folder, background_sub, brightness_thresh);
        u{i,j} = loadInputPattern(data_folder, pattern_blurring);
        %assert( all(u==inputs{i}.Values,'all') )
        
        % get distribution wrt light intensity
        [density_by_input_exp(i,j,:), bins, norm_slope(i,j), c_coeff(i,j), coefficents(i,j,:), agents_by_input(i,j,:), pixels_by_input(i,j,:)] = agentsDensityByInput(inputs{i}.Points, inputs{i}.Values, mask{i,j}, window);
        
        % evaluate quality of fit
        tvd(i,j) = 0.5 * norm(density_by_input_sim(i,:)-squeeze(density_by_input_exp(i,j,:)),1); % Total Variation Distance
        
        % combination of masks over the replicates
        if j==1
            combo_mask{i} = mask{i,1};
        else
            combo_mask{i} = (combo_mask{i}+mask{i,j})>=1;
        end
    end
    
    % weighted average light distribution over the replicates
    mean_dist(i,:) = sum(squeeze(agents_by_input(i,:,:))./squeeze(pixels_by_input(i,:,:)),1);
    mean_dist(i,:) = mean_dist(i,:)/sum(mean_dist(i,:));
end

%% PRINT RESULTS
metrics_of_interest = {tvd}; metrics_tags = ["TVD"]; metrics_color = ['b'];
cmap = linspace2([1,1,1], [1,0.5,0.5], 100)';
x_vec = linspace(window(1),window(2),size(mask,2));
y_vec = linspace(window(3),window(4),size(mask,1));


% multi-exp comparison
% main_fig = figure('Position',[100 100 1900 1000]);
% for i = 1:size(experiments_names,1)  % for each experiment
%     subplot(size(experiments_names,2)+3,size(experiments_names,1),i)
%     box on
%     hold on
%     plotEnvField(inputs{i}.Points, inputs{i}.Values, arena)
%     %plotSwarmInit(xFinal_inWindow{i}, [], inf, inf, Simulation.arena);
%     plotSwarm(xFinal_inWindow{i});
%     axis('equal')
%     axis(window)
%     xticks([])
%     yticks([])
%     title(sim_names{i},'Interpreter','none','FontSize',14)
%     if i==1
%         ylabel('Simulation')
%     end
%     
%     for j = 1:size(experiments_names,2)  % for each replicates
%         x_vec = linspace(window(1),window(2),size(mask{i,j},2));
%         y_vec = linspace(window(3),window(4),size(mask{i,j},1));
%         subplot(size(experiments_names,2)+3,size(experiments_names,1),j*size(experiments_names,1)+i)
%         box on
%         hold on
%         colormap(cmap)
%         imagesc(x_vec,y_vec,u{i,j}')
%         I=imagesc(x_vec,y_vec,cat(3,zeros(size(mask{i,j})),zeros(size(mask{i,j})),mask{i,j}));
%         set(I, 'AlphaData', mask{i,j});
%         axis('equal')
%         axis(window)
%         xticks([])
%         yticks([])
%         title(experiments_names{i},'Interpreter','none')
%         if i==1 && j==1
%             ylabel('Experiments')
%         end
%     end
%     
%     x_vec = linspace(window(1),window(2),size(combo_mask{i},2));
%     y_vec = linspace(window(3),window(4),size(combo_mask{i},1));
%     subplot(size(experiments_names,2)+3,size(experiments_names,1),(j+1)*size(experiments_names,1)+i)
%     box on
%     hold on
%     colormap(cmap)
%     imagesc(x_vec,y_vec,u{i,j}')
%     I=imagesc(x_vec,y_vec,cat(3,zeros(size(mask{i,j})),zeros(size(combo_mask{i})),combo_mask{i}));
%     set(I, 'AlphaData', combo_mask{i});
%     axis('equal')
%     axis(window)
%     xticks([])
%     yticks([])
%     title(experiments_names{i},'Interpreter','none')
%     if i==1
%         ylabel('Combo Experiment')
%     end
%     
%     subplot(size(experiments_names,2)+3,size(experiments_names,1),(j+2)*size(experiments_names,1)+i)
%     hold on
%     b_exp_mean = bar((bins(1:end-1)+bins(2:end))/2,mean_dist(i,:), 1, FaceColor = 'b', FaceAlpha = 0.5);
%     b_sim = bar((bins(1:end-1)+bins(2:end))/2,density_by_input_sim(i,:), 1, FaceColor = 'k', FaceAlpha = 0.4);
%     %[f,xi] = ksdensity(u_values_exp, support=[-0.001,1.001], BoundaryCorrection='reflection');
%     %f=f/sum(f);
%     %plot(xi,f)
%     legend({'REAL','SIMULATED'},'FontSize',14)
%     xlabel('Input intensity','FontSize',14)
%     ylabel('Density of agents','FontSize',14)
%     yticks([0:0.25:1]);
%     text(0.1,max(density_by_input_exp(i,:))*1.10,['TVD=',num2str(tvd(i,j),'%.2f')],'HorizontalAlignment','center','FontSize',14)
%     ylim([0,max(density_by_input_exp(i,:))*1.15])
%     xlim([-0.1,1.1])
%     xticks(round(bins,2))
%     box
%     
% end
% saveas(gcf,fullfile(output_folder, 'multi_exp_comparison_overview'))
% saveas(gcf,fullfile(output_folder, 'multi_exp_comparison_overview'),'png')

% multi-exp comparison
main_fig = figure('Position',[100 100 1900 1000]);
for i = 1:size(experiments_names,1)  % for each experiment
    subplot(4,size(experiments_names,1),i)
    box on
    hold on
    plotEnvField(inputs{i}.Points, inputs{i}.Values, arena)
    %plotSwarmInit(xFinal_inWindow{i}, [], inf, inf, Simulation.arena);
    plotSwarm(xFinal_inWindow{i});
    axis('equal')
    axis(window)
    xticks([])
    yticks([])
    title(sim_names{i},'Interpreter','none','FontSize',14)
    if i==1
        ylabel('Simulation')
    end
    
    x_vec = linspace(window(1),window(2),size(combo_mask{i},2));
    y_vec = linspace(window(3),window(4),size(combo_mask{i},1));
    subplot(4,size(experiments_names,1),size(experiments_names,1)+i)
    box on
    hold on
    colormap(cmap)
    imagesc(x_vec,y_vec,u{i,j}')
    I=imagesc(x_vec,y_vec,cat(3,zeros(size(mask{i,j})),zeros(size(combo_mask{i})),combo_mask{i}));
    set(I, 'AlphaData', combo_mask{i});
    axis('equal')
    axis(window)
    xticks([])
    yticks([])
    title(experiments_names{i},'Interpreter','none')
    if i==1
        ylabel('Combo Experiment')
    end
    
    subplot(4,size(experiments_names,1),2*size(experiments_names,1)+i)
    hold on
    b_exp_mean = bar((bins(1:end-1)+bins(2:end))/2,mean_dist(i,:), 1, FaceColor = 'b', FaceAlpha = 0.5);
    b_sim = bar((bins(1:end-1)+bins(2:end))/2,density_by_input_sim(i,:), 1, FaceColor = 'k', FaceAlpha = 0.4);
    %[f,xi] = ksdensity(u_values_exp, support=[-0.001,1.001], BoundaryCorrection='reflection');
    %f=f/sum(f);
    %plot(xi,f)
    legend({'REAL','SIMULATED'},'FontSize',14)
    xlabel('Input intensity','FontSize',14)
    ylabel('Density of agents','FontSize',14)
    yticks([0:0.25:1]);
    text(0.1,max(density_by_input_exp(i,:))*1.10,['TVD=',num2str(tvd(i,j),'%.2f')],'HorizontalAlignment','center','FontSize',14)
    ylim([0,max(density_by_input_exp(i,:))*1.15])
    xlim([-0.1,1.1])
    xticks(round(bins,2))
    box
    
end
subplot(4,size(experiments_names,1),[1+3*size(experiments_names,1),i+3*size(experiments_names,1)])
hold on
for k=1:length(metrics_of_interest)
    x_pos = [[1:length(tags)]-(length(metrics_of_interest)-1)*0.1+(k-1)*0.2]'-linspace(-1,1,size(experiments_names,2)-1)*0.025;
    plots(k,:)=bar(mean(x_pos,2),metrics_of_interest{k}(:,1),0.15,metrics_color(k),'FaceAlpha',0.5);
    scatter(x_pos,metrics_of_interest{k}(:,2:end),100,metrics_color(k),'MarkerFaceColor','w','LineWidth',1);
    %plots(k,:)=scatter([1:length(tags)]-(length(metrics_of_interest)-1)*0.1+(k-1)*0.2,metrics_of_interest{k}(:,1),100,metrics_color(k),"filled");
end
xticks([1:length(tags)])
xticklabels(tags)
set(gca, 'TickLabelInterpreter', 'none');
xlim([0.7,length(tags)+0.3])
ylim([0, max([metrics_of_interest{:}],[],'all')*1.1])
legend(plots(:,1),metrics_tags,'FontSize',14,'Orientation','horizontal')
box on
set(gca,'XGrid','off','YGrid','on')
saveas(gcf,fullfile(output_folder, 'multi_exp_comparison'))
saveas(gcf,fullfile(output_folder, 'multi_exp_comparison'),'png')


