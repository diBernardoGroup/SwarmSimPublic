clear
close all

defaultParamMicroorg

simulations_folder = '/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations';
simulations_folder = fullfile(simulations_folder,'2024_06_17_GB_absw_noalpha_narrow');
experiments_folder = "/Volumes/DOMEPEN/Experiments";


% half_half
tags = ["half_half"];
sim_names = ["2024_05_30_half_half_1"];
experiments_names = ["2023_06_12_Euglena_2","2023_06_14_Euglena_6","2023_06_15_Euglena_12"];%,"2023_06_26_Euglena_29","2023_06_26_Euglena_30","2023_06_23_Euglena_1","2023_06_23_Euglena_2","2023_06_26_Euglena_2","2023_06_26_Euglena_1"];

% spatial
tags = ["half_half","grad_centr_light","grad_centr_dark","grad_lateral","circle_light","circle_dark"];
sim_names = ["2024_06_06_half_half_1";"2024_06_06_grad_centr_light_1";"2024_06_06_grad_centr_dark_1";"2024_06_06_grad_lateral_1";"2024_06_06_circle_light_1";"2024_06_06_circle_dark_1"];
sim_names = ["2024_06_19_half_half_1";"2024_06_19_grad_centr_light_1";"2024_06_19_grad_centr_dark_1";"2024_06_19_grad_lateral_1";"2024_06_19_circle_light_1";"2024_06_19_circle_dark_1"];
experiments_names = {["2023_06_12_E_2", "2023_06_14_E_6", "2023_06_15_E_12","2023_06_26_E_29","2023_06_26_E_30","2023_06_23_E_1", "2023_06_23_E_2", "2023_06_26_E_2", "2023_06_26_E_1"];
                     ["2023_06_12_E_3", "2023_06_12_E_4", "2023_06_14_E_7", "2023_06_15_E_14","2023_06_23_E_5", "2023_06_23_E_6", "2023_06_26_E_5", "2023_06_26_E_6", "2023_06_26_E_33"];
                     ["2023_06_14_E_10","2023_06_15_E_15","2023_06_23_E_7", "2023_06_23_E_8", "2023_06_23_E_9",  "2023_06_26_E_7","2023_06_26_E_8", "2023_06_26_E_34","2023_06_26_E_35","2023_07_10_E_23","2023_07_10_E_24"];
                     ["2023_06_12_E_5", "2023_06_13_E_16","2023_06_14_E_8", "2023_06_15_E_13","2023_06_23_E_3", "2023_06_23_E_4", "2023_06_26_E_3", "2023_06_26_E_4", "2023_06_26_E_31","2023_06_26_E_32"];
                     ["2023_06_12_E_1", "2023_06_14_E_1", "2023_06_15_E_16","2023_06_23_E_10","2023_06_23_E_11","2023_06_26_E_9", "2023_06_26_E_10","2023_06_26_E_36","2023_06_26_E_37","2023_07_10_E_26"];
                     ["2023_06_13_E_6", "2023_06_13_E_15","2023_06_15_E_17","2023_06_15_E_18","2023_06_23_E_12","2023_06_23_E_13","2023_06_26_E_11","2023_06_26_E_12","2023_06_26_E_38","2023_06_26_E_39","2023_07_10_E_25","2023_07_10_E_22"]};
output_folder = simulations_folder;


%% LOAD DATA
combo_mask = cell(1,length(experiments_names));
assert((length(tags)==length(sim_names)) && (length(sim_names)==length(experiments_names)))
for i = 1:length(experiments_names)  % for each experiment
    
    % load simulation data
    sim_folder = fullfile(simulations_folder,sim_names(i));
    sim_data = load(fullfile(sim_folder,'data.mat'));
    xFinal_inWindow{i} = sim_data.xFinal_inWindow;
    arena = sim_data.Simulation.arena;
    window = [-arena(1),arena(1),-arena(2),arena(2)]/2;
    inputs{i} = sim_data.Environment.Inputs;
    [density_by_input_sim{i}, bins, norm_slope_sim(i), c_coeff_sim(i)] = agentsDensityByInput(inputs{i}.Points, inputs{i}.Values, xFinal_inWindow{i}, window);
    
    % load experiment data
    for j=1:length(experiments_names{i})   % for each replicate
        experiments_names{i}(j) = strrep(experiments_names{i}(j),'_E_','_Euglena_');
        data_folder =  fullfile(experiments_folder,experiments_names{i}(j));
        mask{i}{j} = detectObjects(data_folder, background_sub, brightness_thresh);
        u{i}{j} = loadInputPattern(data_folder, pattern_blurring);
        assert( mean(abs(u{i}{1} - imresize(u{i}{j},size(u{i}{1}))),'all')<0.1, 'Replicates have different inputs' )
        %fprintf('exp %d rep %d size(u)=[%d, %d] mean(u-u1)=%.2f\n',i,j,size(u{i,j}), mean(abs(u{i,1} - imresize(u{i,j},size(u{i,1}))),'all'))
        
        
        % get distribution wrt light intensity
        [density_by_input_exp{i}(j,:), bins, norm_slope{i}(j), c_coeff{i}(j), coefficents{i}(j,:), agents_by_input{i}(j,:), pixels_by_input{i}(j,:)] = agentsDensityByInput(inputs{i}.Points, inputs{i}.Values, mask{i}{j}, window);
        
        % evaluate quality of fit
        tvd{i}(j) = 0.5 * norm(density_by_input_sim{i}-squeeze(density_by_input_exp{i}(j,:)),1); % Total Variation Distance
        
        % combination of masks over the replicates
        if j==1
            combo_mask{i} = mask{i}{1};
        else
            combo_mask{i} = (combo_mask{i}+mask{i}{j})>=1;
        end
    end
    
    % weighted average light distribution over the replicates
    mean_dist{i} = sum(squeeze(agents_by_input{i})./squeeze(pixels_by_input{i}),1);
    mean_dist{i} = mean_dist{i}/sum(mean_dist{i});
    mean_tvd(i) = 0.5 * norm(density_by_input_sim{i}-mean_dist{i},1); % Total Variation Distance
end

%% PRINT RESULTS
metrics_of_interest = {tvd}; metrics_tags = ["TVD"]; metrics_color = ['b'];
cmap = linspace2([1,1,1], [1,0.5,0.5], 100)';
x_vec = linspace(window(1),window(2),size(mask,2));
y_vec = linspace(window(3),window(4),size(mask,1));


% % multi-exp comparison
% main_fig = figure('Position',[100 100 1900 1000]);
% for i = 1:length(experiments_names)  % for each experiment
%     subplot(length(experiments_names{i})+3,length(experiments_names),i)
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
%     for j = 1:length(experiments_names{i})  % for each replicates
%         x_vec = linspace(window(1),window(2),size(mask{i,j},2));
%         y_vec = linspace(window(3),window(4),size(mask{i,j},1));
%         subplot(length(experiments_names{i})+3,length(experiments_names),j*length(experiments_names)+i)
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
%     subplot(length(experiments_names{i})+3,length(experiments_names),(j+1)*length(experiments_names)+i)
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
%     subplot(length(experiments_names{i})+3,length(experiments_names),(j+2)*length(experiments_names)+i)
%     hold on
%     b_exp_mean = bar((bins(1:end-1)+bins(2:end))/2,mean_dist{i}, 1, FaceColor = 'b', FaceAlpha = 0.5);
%     b_sim = bar((bins(1:end-1)+bins(2:end))/2,density_by_input_sim{i}, 1, FaceColor = 'k', FaceAlpha = 0.4);
%     %[f,xi] = ksdensity(u_values_exp, support=[-0.001,1.001], BoundaryCorrection='reflection');
%     %f=f/sum(f);
%     %plot(xi,f)
%     legend({'REAL','SIMULATED'},'FontSize',14)
%     xlabel('Input intensity','FontSize',14)
%     ylabel('Density of agents','FontSize',14)
%     yticks([0:0.25:1]);
%     text(0.1,max(density_by_input_exp{i})*1.10,['TVD=',num2str(mean_tvd(i),'%.2f')],'HorizontalAlignment','center','FontSize',14)
%     ylim([0,max(density_by_input_exp{i})*1.15])
%     xlim([-0.1,1.1])
%     xticks(round(bins,2))
%     box
%     
% end
% saveas(gcf,fullfile(output_folder, 'multi_exp_comparison_overview'))
% saveas(gcf,fullfile(output_folder, 'multi_exp_comparison_overview'),'png')

% multi-exp comparison
main_fig = figure('Position',[100 100 1900 1000]);
for i = 1:length(experiments_names)  % for each experiment
    % simulation final positions
    subplot(4,length(experiments_names),i)
    box on
    hold on
    plotEnvField(inputs{i}.Points, inputs{i}.Values, arena)
    plotSwarm(xFinal_inWindow{i},[],0,inf,inf,false, [], false, 5);
    axis('equal')
    axis(window)
    xticks([])
    yticks([])
    title(sim_names{i},'Interpreter','none','FontSize',12)
    if i==1
        ylabel('Simulation','FontSize',12)
    end
    
    % combo experiment mask
    x_vec = linspace(window(1),window(2),size(combo_mask{i},2));
    y_vec = linspace(window(3),window(4),size(combo_mask{i},1));
    subplot(4,length(experiments_names),length(experiments_names)+i)
    box on
    hold on
    colormap(cmap)
    imagesc(x_vec,y_vec,u{i}{1}')
    I=imagesc(x_vec,y_vec,cat(3,zeros(size(mask{i}{1})),zeros(size(combo_mask{i})),combo_mask{i}));
    set(I, 'AlphaData', combo_mask{i});
    axis('equal')
    axis(window)
    xticks([])
    yticks([])
    if i==1
        ylabel('Combo Experiment','FontSize',12)
    end
    
    % distributions wrt light
    subplot(4,length(experiments_names),2*length(experiments_names)+i)
    hold on
    b_exp_mean = bar((bins(1:end-1)+bins(2:end))/2,mean_dist{i}, 1, FaceColor = 'b', FaceAlpha = 0.5);
    b_sim = bar((bins(1:end-1)+bins(2:end))/2,density_by_input_sim{i}, 1, FaceColor = 'k', FaceAlpha = 0.4);
    %[f,xi] = ksdensity(u_values_exp, support=[-0.001,1.001], BoundaryCorrection='reflection');
    %f=f/sum(f);
    %plot(xi,f)
    if i==length(experiments_names)
    legend({'REAL','SIMULATED'},'FontSize',12)%,'Location','best')
    end
    xlabel('Input intensity','FontSize',12)
    ylabel('Density of agents','FontSize',12)
    yticks([0:0.25:1]);
    text(0,0.85,['TVD=',num2str(mean_tvd(i),'%.2f')],'FontSize',12)%,'HorizontalAlignment','center'
    ylim([0,1])
    xlim([-0.1,1.1])
    xticks(round(bins,2))
    box
    
end
% metrics
subplot(4,length(experiments_names),[1+3*length(experiments_names),i+3*length(experiments_names)])
hold on
for k=1:length(metrics_of_interest)
    for i=1:length(experiments_names)
        %x_pos = [[1:length(tags)]-(length(metrics_of_interest)-1)*0.1+(k-1)*0.2]'+linspace(-1,1,length(experiments_names{i}))*0.05;
        x_pos = [i-(length(metrics_of_interest)-1)*0.1+(k-1)*0.2]'+linspace(-1,1,length(experiments_names{i}))*0.05;
        plots(k,i)=bar(mean(x_pos),mean_tvd(i),0.15,metrics_color(k),'FaceAlpha',0.5);
        scatter(x_pos, metrics_of_interest{k}{i},100,metrics_color(k),'MarkerFaceColor','w','LineWidth',1.25);
    end
    %plots(k,:)=scatter([1:length(tags)]-(length(metrics_of_interest)-1)*0.1+(k-1)*0.2,metrics_of_interest{k}(:,1),100,metrics_color(k),"filled");
    %text(length(tags)+0.4, max([metrics_of_interest{:}],[],'all')*(0.7-k*0.1), ['med ',char(metrics_tags(k)),'=',num2str(median(mean_tvd),'%.2f')],'FontSize',12)
end
xticks([1:length(tags)])
xticklabels(tags)
set(gca, 'TickLabelInterpreter', 'none');
set(gca,'FontSize',12)
xlim([0.7,length(tags)+0.3])
ylim([0, 0.5])
legend(plots(:,1),metrics_tags,'FontSize',12,'Orientation','horizontal')
box on
set(gca,'XGrid','off','YGrid','on')
saveas(gcf,fullfile(output_folder, 'multi_exp_comparison_spatial'))
saveas(gcf,fullfile(output_folder, 'multi_exp_comparison_spatial'),'png')


