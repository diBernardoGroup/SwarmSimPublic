%
% CombineResults
%
%   See also: CombineResultsSpatial.
%
%   Authors:    Andrea Giusti and Davide Salzano
%   Date:       2024
%

clear
close all

defaultParamMicroorg;

% ANDREA
simulations_folder = 'Output';
experiments_folder = "/PATH/TO/THE/EXPERIMENTS";         % path to the DOME experiments data


%% EUGLENA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

simulations_folder = fullfile(simulations_folder,'2024_06_17_E_GB_absw_noalpha_narrow','2024_08_02_time_2k_v2'); % select subfolder

% % switch 10s all
% tags = "switch_10";
% sim_names = "2024_04_17_switch_10_2";
% sim_names = "2024_05_09_switch_10_3";
% experiments_names = ["comparisons/Euglena_switch_10/combo5", "2023_06_15_Euglena_7", "2023_06_26_Euglena_23","2023_06_26_Euglena_24", "2023_07_10_Euglena_15", "2023_07_10_Euglena_16"];
%
% % switch 10s selected
% tags = "switch_10";
% sim_names = "2024_04_17_switch_10_2";
% sim_names = "2024_05_09_switch_10_3";
% experiments_names = ["comparisons/Euglena_switch_10/combo3","2023_06_15_Euglena_7", "2023_06_26_Euglena_24", "2023_07_10_Euglena_15"];
%
% % switch 5s
% tags = "switch_5";
% sim_names = "2024_04_17_switch_5_1";
% sim_names = "2024_05_09_switch_5_1";
% experiments_names = ["comparisons/Euglena_switch_5/combo", "2023_06_15_Euglena_8", "2023_06_26_Euglena_25", "2023_07_10_Euglena_18"];
%
% % switch 1s
% tags = "switch_1";
% sim_names = "2024_04_17_switch_1_2";
% sim_names = "2024_05_09_switch_1_1";
% experiments_names = ["comparisons/Euglena_switch_1/combo", "2023_06_15_Euglena_11", "2023_06_26_Euglena_28", "2023_07_10_Euglena_19"];
%
% % off
% tags = "OFF";
% sim_names = "2024_04_17_OFF_2";
% sim_names = "2024_05_09_OFF_2";
% experiments_names = ["comparisons/Euglena_OFF/combo", "2023_06_15_Euglena_1", "2023_06_26_Euglena_13", "2023_07_10_Euglena_6"];
%
% OFF-ON-OFF 75
% tags = "75_ON";
% sim_names = "2024_05_09_75_ON_2";
% experiments_names = ["comparisons/Euglena_75_ON/combo_old", "2023_06_15_Euglena_2", "2023_06_26_Euglena_16", "2023_07_10_Euglena_8"]; %old combo
% experiments_names = ["comparisons/Euglena_75_ON/combo", "2023_06_15_Euglena_2", "2023_06_26_Euglena_15", "2023_07_10_Euglena_8"];
%
% % OFF-ON-OFF 150
% tags = "150_ON";
% sim_names = "2024_04_17_150_ON_3";
% sim_names = "2024_05_09_150_ON_2";
% experiments_names = ["comparisons/Euglena_150_ON/combo", "2023_06_15_Euglena_3", "2023_06_26_Euglena_18", "2023_07_10_Euglena_10"];
%
% % OFF-ON-OFF 255
% tags = "255_ON";
% sim_names = "2024_04_17_255_ON_3";
% sim_names = "2024_05_09_255_ON_3";
% experiments_names = ["comparisons/Euglena_255_ON/combo", "2023_06_15_Euglena_4", "2023_06_26_Euglena_20", "2023_07_10_Euglena_12"];
%
% % Ramp
% tags = "ramp";
% sim_names = "2024_04_17_ramp_2";
% sim_names = "2024_05_09_ramp_1";
% experiments_names = ["comparisons/Euglena_ramp/combo", "2023_06_15_Euglena_6", "2023_06_26_Euglena_22", "2023_07_10_Euglena_14"];

% all
tags = ["switch_10","switch_5","switch_1","ramp", "OFF","75_ON", "150_ON", "255_ON"];
sim_names = ["experiment_switch_10_1", "experiment_switch_5_1", "experiment_switch_1_1", "experiment_ramp_1", "experiment_OFF_1", "experiment_75_ON_1", "experiment_150_ON_1", "experiment_255_ON_1"]';
experiments_names = [
    "comparisons/Euglena_switch_10/combo3","2023_06_15_Euglena_7", "2023_06_26_Euglena_24", "2023_07_10_Euglena_15";
    "comparisons/Euglena_switch_5/combo", "2023_06_15_Euglena_8", "2023_06_26_Euglena_25", "2023_07_10_Euglena_18";
    "comparisons/Euglena_switch_1/combo", "2023_06_15_Euglena_11", "2023_06_26_Euglena_28", "2023_07_10_Euglena_19";
    "comparisons/Euglena_ramp/combo", "2023_06_15_Euglena_6", "2023_06_26_Euglena_22", "2023_07_10_Euglena_14";
    "comparisons/Euglena_OFF/combo", "2023_06_15_Euglena_1", "2023_06_26_Euglena_13", "2023_07_10_Euglena_6";
    "comparisons/Euglena_75_ON/combo", "2023_06_15_Euglena_2", "2023_06_26_Euglena_15", "2023_07_10_Euglena_8";
    "comparisons/Euglena_150_ON/combo", "2023_06_15_Euglena_3", "2023_06_26_Euglena_18", "2023_07_10_Euglena_10";
    "comparisons/Euglena_255_ON/combo", "2023_06_15_Euglena_4", "2023_06_26_Euglena_20", "2023_07_10_Euglena_12"];
output_folder = simulations_folder;

% scenario duration
% tags = tags(1:4);
% sim_names = sim_names(1:4);
% experiments_names = experiments_names(1:4,:);

% scenario intensity
% tags = tags(5:8);
% sim_names = sim_names(5:8);
% experiments_names = experiments_names(5:8,:);


%% VOLVOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % simulations_folder = fullfile(simulations_folder,'2024_08_01_V_GB_meaninit','2024_08_01_time_500');% select subfolder
% simulations_folder = fullfile(simulations_folder,'2024_08_02_V_GB_60s','2024_08_02_time_500');     % select subfolder
% 
% % all
% tags = ["switch_10","switch_5","switch_1","ramp", "OFF","75_ON", "150_ON", "255_ON"];
% sim_names = ["experiment_switch_10_1", "experiment_switch_5_1", "experiment_switch_1_1", "experiment_ramp_1", "experiment_OFF_1", "experiment_75_ON_1", "experiment_150_ON_1", "experiment_255_ON_1"]';
% experiments_names = ["comparisons/Volvox_switch_10/combo3","2023_07_04_V_13","2023_07_05_V_2","2023_07_06_V_5";
%     "comparisons/Volvox_switch_5/combo", "2023_07_04_V_14","2023_07_05_V_8","2023_07_06_V_12";
%     "comparisons/Volvox_switch_1/combo", "2023_07_04_V_16","2023_07_05_V_9","2023_07_06_V_13";
%     "comparisons/Volvox_ramp/combo", "2023_07_04_V_11", "2023_07_05_V_6", "2023_07_06_V_14";
%     "comparisons/Volvox_OFF/combo", "2023_07_04_V_2", "2023_07_05_V_1", "2023_07_06_V_10";
%     "comparisons/Volvox_75_ON/combo", "2023_07_04_V_4", "2023_07_05_V_3", "2023_07_06_V_18";
%     "comparisons/Volvox_150_ON/combo", "2023_07_04_V_7", "2023_07_05_V_4", "2023_07_06_V_20";
%     "comparisons/Volvox_255_ON/combo", "2023_07_04_V_9", "2023_07_05_V_5", "2023_07_06_V_23"];
% output_folder = simulations_folder;

deltaT = 0.5;
timeInstants = [0:deltaT:180];


%% LOAD DATA

assert((length(tags)==length(sim_names)) && (length(sim_names)==size(experiments_names,1)))
for i = 1:size(experiments_names,1)  % for each experiment
    
    % load simulation data
    sim_folder = fullfile(simulations_folder,sim_names(i));
    sim_data = load(fullfile(sim_folder,'data.mat'));
    speed_sim{i} = sim_data.speed;
    omega_sim{i} = sim_data.omega;
    u{i} = sim_data.u;
        
    % load experiment data
    for j=1:size(experiments_names,2)   % for each replicate
        experiments_names(i,j) = strrep(experiments_names(i,j),'_E_','_Euglena_');
        experiments_names(i,j) = strrep(experiments_names(i,j),'_V_','_Volvox_');
        if j==1
            data_folder =  fullfile(experiments_folder,experiments_names(i,j));
        else
            subFolfder = getLastTracking(fullfile(experiments_folder,experiments_names(i,j)));
            data_folder =  fullfile(experiments_folder,experiments_names(i,j),subFolfder);
        end
        speeds{i,j}  = load(fullfile(data_folder,'speeds_smooth.txt'));
        omegas{i,j}  = load(fullfile(data_folder,'ang_vel_smooth.txt'));
        inputs = load(fullfile(data_folder,'inputs.txt'));
        assert( all(inputs(:,1)/255==u{i},'all'), 'Experiment has inputs different from simulation!' )
        
        % evaluate quality of fit
        [metrics_of_interest] = compareResults({data_folder,sim_folder}, sim_folder, false, Render);
        wmape_speed(i,j) = metrics_of_interest{1};
        wmape_omega(i,j) = metrics_of_interest{2};
        wmape_total(i,j) = metrics_of_interest{3};
    end
    
    % evaluate quality of fit
%     for j=1:size(experiments_names,2)   % for each replicate
%         overlap = min(size(omega_sim{i},1),size(omegas{i,j},1));
%         nmse_speed(i,j) = goodnessOfFit(median(speed_sim{i},2,'omitnan'), median(speeds{i,j},2,'omitnan'), 'NMSE');
%         nmse_omega(i,j) = goodnessOfFit(median(abs(omega_sim{i}(1:overlap,:)),2,'omitnan'), median(abs(omegas{i,j}(1:overlap,:)),2,'omitnan'), 'NMSE');
%         nmse_total(i,j) = mean([nmse_speed(i,j), nmse_omega(i,j)]);
%      
%         mape_speed(i,j) = mape(median(speed_sim{i},2,'omitnan'), median(speeds{i,j},2,'omitnan'));
%         mape_omega(i,j) = mape(median(abs(omega_sim{i}(1:overlap,:)),2,'omitnan'), median(abs(omegas{i,j}(1:overlap,:)),2,'omitnan'));
%         mape_total(i,j) = mean([mape_speed(i,j), mape_omega(i,j)]);
%     
%         wmape_speed(i,j) = mape(median(speed_sim{i},2,'omitnan'), median(speeds{i,j},2,'omitnan'),'wMAPE');
%         wmape_omega(i,j) = mape(median(abs(omega_sim{i}(1:overlap,:)),2,'omitnan'), median(abs(omegas{i,j}(1:overlap,:)),2,'omitnan'),'wMAPE');
%         wmape_total(i,j) = mean([wmape_speed(i,j), wmape_omega(i,j)]);
%     end
end

%% PRINT RESULTS
metrics_of_interest = {wmape_speed, wmape_total, wmape_omega}; metrics_tags = ["wmape_v", "wmape_{tot}", "wmape\omega"]; metrics_color = ['b','k','r'];


%%%%%%%%%%%%%%%%%%%%%%%%% Single experiment plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% for i = 1:size(experiments_names,1)  % for each experiment
% if size(experiments_names,1) == 1
%     fileID = fopen(fullfile(simulations_folder, sim_names(i), 'multi_exp_comparison.txt'),'wt');
%     fprintf(fileID,'DOME multi experiment comparison\n\n');
%     fprintf(fileID,'Date: %s\n',datestr(now, 'dd/mm/yy'));
%     fprintf(fileID,'Time: %s\n\n',datestr(now, 'HH:MM'));
%     fprintf(fileID,'Experiment\t\tNMSE speed\tNMSE omega\tNMSE tot\n');
%     for j=1:size(experiments_names,2)
%         fprintf(fileID,'%s\t',experiments_names(i,j));
%         fprintf(fileID,'%.2f\t\t%.2f\t\t%.2f\n',nmse_speed(i,j),nmse_omega(i,j),nmse_total(i,j));
%     end
%     fclose(fileID);
%
%     figure %time plot
%     subplot(2,1,1)
%     hold on
%     ylim([0,100])
%     highlightInputs(timeInstants, u{i}, 'r', 0.25)
%     for j=2:size(experiments_names,2)
%         plot(timeInstants, median(speeds{i,j},2,'omitnan'),'b', color=[0.5,0.5,1]);
%     end
%     l1=plot(timeInstants, median(speeds{i,1},2,'omitnan'),'b',LineWidth=2);
%     l2=plot(timeInstants, median(speed_sim{i},2,'omitnan'),'k',LineWidth=2);
%     xlim([0,max(timeInstants)])
%     xlabel('$t$ [s]','Interpreter','Latex','FontSize',16)
%     ylabel('$v$ [px/s]','Interpreter','Latex','FontSize',16)
%     legend([l1,l2],'REAL','SIMULATED')
%     box on
%     subplot(2,1,2)
%     hold on
%     ylim([0,1.5])
%     highlightInputs(timeInstants, u{i}, 'r', 0.25)
%     for j=2:size(experiments_names,2)
%         plot(timeInstants, median(abs(omegas{i,j}),2,'omitnan'),'b', color=[0.5,0.5,1]);
%     end
%     l1=plot(timeInstants, median(abs(omegas{i,1}),2,'omitnan'),'b',LineWidth=2);
%     l2=plot(timeInstants, median(abs(omega_sim{i}),2,'omitnan'),'k',LineWidth=2);
%     xlabel('$t$ [s]','Interpreter','Latex','FontSize',16)
%     ylabel('$|\omega|$ [rad/s]','Interpreter','Latex','FontSize',16)
%     legend([l1,l2],'REAL','SIMULATED')
%     box on
%     saveas(gcf,fullfile(simulations_folder,sim_names(i), 'multi_exp_comparison_time_plot'))
%     saveas(gcf,fullfile(simulations_folder,sim_names(i), 'multi_exp_comparison_time_plot'),'png')
%
%
%     figure % metrics scatter single exp
%     hold on
%     for k=1:length(metrics_of_interest)
%         x_pos = [1-(length(metrics_of_interest)-1)*0.1+(k-1)*0.2]'-linspace(-1,1,size(experiments_names,2)-1)*0.033;
%         plots(k,:)=bar(mean(x_pos,2),metrics_of_interest{k}(:,1),0.15,metrics_color(k),'FaceAlpha',0.5);
%         scatter(x_pos,metrics_of_interest{k}(:,2:end),100,metrics_color(k),'MarkerFaceColor','w','LineWidth',1);
%     end
%     xticks([1])
%     xticklabels(tags(i))
%     set(gca, 'TickLabelInterpreter', 'none');
%     xlim([0.5,1.5])
%     all_metrics = [metrics_of_interest{:}];
%     ylim([0, max(all_metrics(i,:),[],'all')*1.1])
%     legend(plots(:,1),metrics_tags)
%     set(gca,'FontSize',14)
%     box on
%     set(gca,'XGrid','off','YGrid','on')
%     saveas(gcf,fullfile(fullfile(simulations_folder,sim_names(i)), 'multi_exp_comparison_NMSE'))
%     saveas(gcf,fullfile(fullfile(simulations_folder,sim_names(i)), 'multi_exp_comparison_NMSE'),'png')
%
% end


%%%%%%%%%%%%%%%%%%%%%%%% MULTIPLE EXP COMPARISON (subplot) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if size(experiments_names,1) > 1 % multi-exp comparison
    main_fig = figure('Position',[100 100 1900 600]);
    for i = 1:size(experiments_names,1)  % for each experiment
        subplot(3,size(experiments_names,1),i)
        hold on
        ylim([0,100])
        highlightInputs(timeInstants, u{i}, 'r', 0.25)
        for j=2:size(experiments_names,2)
            plot(timeInstants, median(speeds{i,j},2,'omitnan'),'b', color=[0.5,0.5,1]);
        end
        l1=plot(timeInstants, median(speeds{i,1},2,'omitnan'),'b',LineWidth=2);
        l2=plot(timeInstants, median(speed_sim{i},2,'omitnan'),'k',LineWidth=2);
        xlim([0,max(timeInstants)])
        xticks(linspace(0,max(timeInstants),4))
        if i==size(experiments_names,1); legend([l1,l2],'REAL','SIMULATED'); end
        box on
        subplot(3,size(experiments_names,1),i+size(experiments_names,1))
        hold on
        ylim([0,1.5])
        highlightInputs(timeInstants, u{i}, 'r', 0.25)
        for j=2:size(experiments_names,2)
            overlap = min(size(omega_sim{i},1),size(omegas{i,j},1));
            plot(timeInstants, median(abs(omegas{i,j}),2,'omitnan'),'b', color=[0.5,0.5,1]);
        end
        overlap = min(size(omega_sim{i},1),size(omegas{i,1},1));
        l1=plot(timeInstants(1:overlap), median(abs(omegas{i,1}(1:overlap,:)),2,'omitnan'),'b',LineWidth=2);
        l2=plot(timeInstants, median(abs(omega_sim{i}),2,'omitnan'),'k',LineWidth=2);
        xlim([0,max(timeInstants)])
        xticks(linspace(0,max(timeInstants),4))
        xlabel('$t$ [s]','Interpreter','Latex','FontSize',16)
        %ylabel('$|\omega|$ [rad/s]','Interpreter','Latex','FontSize',16)
        if i==size(experiments_names,1); legend([l1,l2],'REAL','SIMULATED'); end
        box on
    end
    
    % bottom panel - metrics
    subplot(3,size(experiments_names,1),1); ylabel('$v$ [$\mu$m/s]','Interpreter','Latex','FontSize',14)
    subplot(3,size(experiments_names,1),1+size(experiments_names,1)); ylabel('$|\omega|$ [rad/s]','Interpreter','Latex','FontSize',14);
    subplot(3,size(experiments_names,1),[1+2*size(experiments_names,1),i+2*size(experiments_names,1)])
    hold on
    for k=1:length(metrics_of_interest)
        x_pos = [[1:length(tags)]-(length(metrics_of_interest)-1)*0.1+(k-1)*0.2]'-linspace(-1,1,size(experiments_names,2))*0.025;
        plots(k,:)=bar(mean(x_pos,2),metrics_of_interest{k}(:,1),0.15,metrics_color(k),'FaceAlpha',0.5);
        scatter(x_pos(:,2:end),metrics_of_interest{k}(:,2:end),100,metrics_color(k),'MarkerFaceColor','w','LineWidth',1.25);
        text(length(tags)+0.4, 1.2*(0.7-k*0.1), ['med ',char(metrics_tags(k)),'=',num2str(median(metrics_of_interest{k}(:,1)),'%.2f')],'FontSize',12)
    end
    xticks([1:length(tags)])
    xticklabels(tags)
    set(gca, 'TickLabelInterpreter', 'none');
    xlim([0.7,length(tags)+0.3])
    ylim([0, max([metrics_of_interest{:}],[],'all')*1.1])
    ylim([0,1.2])
    legend(plots(:,1),metrics_tags,'FontSize',14,'Orientation','horizontal')
    box on
    set(gca,'XGrid','off','YGrid','on')
    saveas(gcf,fullfile(output_folder, 'multi_exp_comparison'))
    saveas(gcf,fullfile(output_folder, 'multi_exp_comparison'),'png')

end

% Metrics only
figure('Position',[100 100 1900 400]);
% metrics_of_interest = metrics_of_interest(2);
% metrics_tags= metrics_tags(2);
hold on
for k=1:length(metrics_of_interest)
    x_pos = [[1:length(tags)]-(length(metrics_of_interest)-1)*0.1+(k-1)*0.2]'-linspace(-1,1,size(experiments_names,2))*0.025;
    plots(k,:)=bar(mean(x_pos,2),metrics_of_interest{k}(:,1),0.15,metrics_color(k),'FaceAlpha',0.5);
    scatter(x_pos(:,2:end),metrics_of_interest{k}(:,2:end),100,metrics_color(k),'MarkerFaceColor','w','LineWidth',1.25);
    text(length(tags)+0.4, 1.2*(0.7-k*0.1), ['med ',char(metrics_tags(k)),'=',num2str(median(metrics_of_interest{k}(:,1)),'%.2f')],'FontSize',12)
end
xticks([1:length(tags)])
xticklabels(tags)
set(gca, 'TickLabelInterpreter', 'none');
xlim([0.7,length(tags)+0.3])
ylim([0, max([metrics_of_interest{:}],[],'all')*1.1])
ylim([0,1.2])
legend(plots(:,1),metrics_tags,'FontSize',14,'Orientation','horizontal')
box on
set(gca,'XGrid','off','YGrid','on')
fig=gcf; fig.Units = fig.PaperUnits; fig.PaperSize = fig.Position(3:4); % set correct pdf size
saveas(gcf,fullfile(output_folder, 'multi_exp_comparison_metrics'))
saveas(gcf,fullfile(output_folder, 'multi_exp_comparison_metrics'),'pdf')

%%%%%%%%%%%%%%%%%%%%%%%% MULTIPLE EXP COMPARISON (FIGURES) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% if size(experiments_names,1) > 1 % multi-exp comparison
%     
%     
%     for i = 1:size(experiments_names,1)  % for each experiment
%         
%         %Defining the position of the image (row and col)
%         win = get(0,'screensize');
%         wid = win(3)/3;
%         hig = 300;
%         col = mod((i-1)*wid,win(3));
%         row = floor((i-1)*wid/win(3));
%         figure('Position',[col row*hig  wid hig]);
%         
%         %SPEED PLOT
%         subplot(2,1,1);
%         ylim([0,100])
%         %Draw the inputs
%         highlightInputs(timeInstants, u{i}, Render.cmap_inputs, 0.7);
%         %Plot the single replicates
%         for j=2:size(experiments_names,2)
%             plot(timeInstants, median(speeds{i,j},2,'omitnan'), color=Render.exp_c);
%         end
%         %Plot the combo and the simulated trajecotry
%         l1=plot(timeInstants, median(speeds{i,1},2,'omitnan'),LineWidth=2,color=Render.exp_c);
%         l2=plot(timeInstants, median(speed_sim{i},2,'omitnan'),LineWidth=2,color=Render.sim_c);
%         %Set the axes limits
%         xlim([0,max(timeInstants)])
%         xticks(linspace(0,max(timeInstants),4))
%         xlabel('$t$ [s]','Interpreter','Latex','FontSize',16)
%         ylabel('$v$ [$\mu$m/s]','Interpreter','Latex','FontSize',16)
%         if i==size(experiments_names,1); legend([l1,l2],'REAL','SIMULATED'); end
%         box on;
%         
%         %ANGUALR VALOCITY
%         subplot(2,1,2);
%         hold on
%         ylim([0,1.5])
%         %Draw the inputs
%         highlightInputs(timeInstants, u{i}, Render.cmap_inputs, 0.7);
%         %Plot the single replicates
%         for j=2:size(experiments_names,2)
%             overlap = min(size(omega_sim{i},1),size(omegas{i,j},1));
%             plot(timeInstants, median(abs(omegas{i,j}),2,'omitnan'), color=Render.exp_c);
%         end
%         %Plot the combo and the simulated trajecotry
%         overlap = min(size(omega_sim{i},1),size(omegas{i,1},1));
%         l1=plot(timeInstants(1:overlap), median(abs(omegas{i,1}(1:overlap,:)),2,'omitnan'),LineWidth=2, color=Render.exp_c);
%         l2=plot(timeInstants, median(abs(omega_sim{i}),2,'omitnan'),LineWidth=2, color=Render.sim_c);
%         %Set the axes limits
%         xlim([0,max(timeInstants)])
%         xticks(linspace(0,max(timeInstants),4))
%         xlabel('$t$ [s]','Interpreter','Latex','FontSize',16)
%         ylabel('$|\omega|$ [rad/s]','Interpreter','Latex','FontSize',16)
%         if i==size(experiments_names,1); legend([l1,l2],'REAL','SIMULATED'); end
%         box on
%         
%     end
% end

