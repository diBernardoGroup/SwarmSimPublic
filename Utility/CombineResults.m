clear
close all

simulations_folder = '/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations';
experiments_folder = "/Volumes/DOMEPEN/Experiments";

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
tags = ["switch_10","switch_5","switch_1","75_ON", "150_ON", "255_ON", "OFF", "ramp"];
%sim_names = ["2024_04_17_switch_10_2", "2024_04_17_switch_5_1", "2024_04_17_switch_1_2", "2024_04_17_75_ON_1", "2024_04_17_150_ON_3", "2024_04_17_255_ON_3", "2024_04_17_OFF_2", "2024_04_17_ramp_2"]';
sim_names = ["2024_05_09_switch_10_3", "2024_05_09_switch_5_1", "2024_05_09_switch_1_1", "2024_05_09_75_ON_2", "2024_05_09_150_ON_2", "2024_05_09_255_ON_3", "2024_05_09_OFF_2", "2024_05_09_ramp_1"]';
sim_names = ["2024_06_03_switch_10_1", "2024_06_03_switch_5_1", "2024_06_03_switch_1_1", "2024_06_03_75_ON_1", "2024_06_03_150_ON_1", "2024_06_03_255_ON_1", "2024_06_03_OFF_1", "2024_06_03_ramp_1"]'; % manual tuning
experiments_names = ["comparisons/Euglena_switch_10/combo3","2023_06_15_Euglena_7", "2023_06_26_Euglena_24", "2023_07_10_Euglena_15";
    "comparisons/Euglena_switch_5/combo", "2023_06_15_Euglena_8", "2023_06_26_Euglena_25", "2023_07_10_Euglena_18";
    "comparisons/Euglena_switch_1/combo", "2023_06_15_Euglena_11", "2023_06_26_Euglena_28", "2023_07_10_Euglena_19";
    "comparisons/Euglena_75_ON/combo", "2023_06_15_Euglena_2", "2023_06_26_Euglena_15", "2023_07_10_Euglena_8";
    "comparisons/Euglena_150_ON/combo", "2023_06_15_Euglena_3", "2023_06_26_Euglena_18", "2023_07_10_Euglena_10";
    "comparisons/Euglena_255_ON/combo", "2023_06_15_Euglena_4", "2023_06_26_Euglena_20", "2023_07_10_Euglena_12";
    "comparisons/Euglena_OFF/combo", "2023_06_15_Euglena_1", "2023_06_26_Euglena_13", "2023_07_10_Euglena_6";
    "comparisons/Euglena_ramp/combo", "2023_06_15_Euglena_6", "2023_06_26_Euglena_22", "2023_07_10_Euglena_14"];
output_folder = '/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/comparison/Euglena input/manual tuning';

% scenario duration
% tags = ["switch_10","switch_5","switch_1","ramp"];
% sim_names = ["2024_05_28_switch_10_1", "2024_05_28_switch_5_1", "2024_05_28_switch_1_1", "2024_05_28_ramp_1"]';
% experiments_names = ["comparisons/Euglena_switch_10/combo3","2023_06_15_Euglena_7", "2023_06_26_Euglena_24", "2023_07_10_Euglena_15";
%                      "comparisons/Euglena_switch_5/combo", "2023_06_15_Euglena_8", "2023_06_26_Euglena_25", "2023_07_10_Euglena_18";
%                      "comparisons/Euglena_switch_1/combo", "2023_06_15_Euglena_11", "2023_06_26_Euglena_28", "2023_07_10_Euglena_19";
%                      "comparisons/Euglena_ramp/combo", "2023_06_15_Euglena_6", "2023_06_26_Euglena_22", "2023_07_10_Euglena_14"];
% output_folder = '/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/comparison/Euglena input/Scenario duration';

% % scenario intensity
% tags = [ "OFF", "75_ON", "150_ON", "255_ON"];
% sim_names = ["2024_05_28_OFF_1", "2024_05_28_75_ON_1", "2024_05_28_150_ON_1", "2024_05_28_255_ON_1"]';
% experiments_names = ["comparisons/Euglena_OFF/combo", "2023_06_15_Euglena_1", "2023_06_26_Euglena_13", "2023_07_10_Euglena_6";
%                     "comparisons/Euglena_75_ON/combo", "2023_06_15_Euglena_2", "2023_06_26_Euglena_15", "2023_07_10_Euglena_8";
%                     "comparisons/Euglena_150_ON/combo", "2023_06_15_Euglena_3", "2023_06_26_Euglena_18", "2023_07_10_Euglena_10";
%                     "comparisons/Euglena_255_ON/combo", "2023_06_15_Euglena_4", "2023_06_26_Euglena_20", "2023_07_10_Euglena_12"];
% output_folder = '/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/comparison/Euglena input/Scenario intensity';
  
% % compare identifications
% tags = [ "oldexp_oldsim_oldid", "newexp_newsim_oldid", "newid_BE_BE_ds1", "newid_BE_BE_ds2", "newid_BE_BE_ds3", "newid_BE_grad_ds1", "newid_BE_delgrad_ds1", "newid_BE_noalpha"];
% sim_names = ["2024_04_17_switch_10_2", "2024_05_08_switch_10_4", "2024_05_08_switch_10_6", "2024_05_08_switch_10_8", "2024_05_08_switch_10_9", "2024_05_08_switch_10_5", "2024_05_13_switch_10_1", "2024_05_13_switch_10_2"]';
% experiments_names = ["comparisons/Euglena_switch_10/combo5_old";"comparisons/Euglena_switch_10/combo5"; "comparisons/Euglena_switch_10/combo5"; "comparisons/Euglena_switch_10/combo5"; "comparisons/Euglena_switch_10/combo5"; "comparisons/Euglena_switch_10/combo5"; "comparisons/Euglena_switch_10/combo5"; "comparisons/Euglena_switch_10/combo5"];
% output_folder = '/Volumes/DOMEPEN/Experiments/comparisons/Euglena_switch_10/combo5';

% % compare identifications
% tags = [ "signed", "nosign", "signed_nomu", "nosign_nomu"];
% sim_names = ["OLS+GB_ds1_diff_sign", "OLS+GB_ds1_diff_nosign", "OLS+GB_ds1_diff_sign_nomu", "OLS+GB_ds1_diff_nosign_nomu"]';
% experiments_names = ["comparisons/Euglena_switch_10/combo5"; "comparisons/Euglena_switch_10/combo5"; "comparisons/Euglena_switch_10/combo5"; "comparisons/Euglena_switch_10/combo5"];
% simulations_folder = '/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/comparison/Identifications';
% output_folder = simulations_folder;

% sim_names =    ["identification_GB_nolim";      "identification_GB_nolim";
%                 "identification_GB_nolim_nomu"; "identification_GB_nolim_nomu";
%                 "identification_GB_lim";        "identification_GB_lim";
%                 "identification_GB_lim_nomu";   "identification_GB_lim_nomu";
%                 ];
% experiments_names = repelem(["comparisons/Euglena_switch_10/combo5"], length(sim_names))';
% tags = ["GA_nolim","GB_nolim","GA_nolim_nomu","GB_nolim_nomu","GA_lim","GB_lim","GA_lim_nomu","GB_lim_nomu"];
% simulations_folder = '/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/comparison/Identifications';
% output_folder = '/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/comparison/Identifications/id_comparison_GA_GB';

% sim_names =    [ "identification_GB_lim"; "identification_GB_lim_nomu";
%                 "identification_GB_lim_v_nomu"; "identification_GB_lim_w_nomu";
%                 "identification_GB_lim_b_nomu";
%                 "identification_GB_lim_bv_nomu"; "identification_GB_lim_bw_nomu";
%                 "identification_GB_lim_b"; "identification_GB_lim_b_discardmu"
%                 ];
% tags = ["GB_lim","GB_lim_nomu","GB_lim_v_nomu","GB_lim_w_nomu","GB_lim_b_nomu","GB_lim_bv_nomu","GB_lim_bw_nomu","GB_lim_b","GB_lim_b_discardmu"];
% experiments_names = repelem(["comparisons/Euglena_switch_10/combo5"], length(sim_names))';
% simulations_folder = '/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/comparison/Identifications';
% output_folder = '/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/comparison/Identifications/id_comparison_GA_GB';

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
        if j==1
            data_folder =  fullfile(experiments_folder,experiments_names(i,j));
        else
            subFolfder = getLastTracking(fullfile(experiments_folder,experiments_names(i,j)));
            data_folder =  fullfile(experiments_folder,experiments_names(i,j),subFolfder);
        end
        speeds{i,j}  = load(fullfile(data_folder,'speeds_smooth.txt'));
        omegas{i,j}  = load(fullfile(data_folder,'ang_vel_smooth.txt'));
        
        inputs = load(fullfile(data_folder,'inputs.txt'));
        assert( all(inputs(:,1)/255==u{i},'all') )
    end
    
    % evaluate quality of fit
    for j=1:size(experiments_names,2)   % for each replicate
        overlap = min(size(omega_sim{i},1),size(omegas{i,j},1));
        nmse_speed(i,j) = goodnessOfFit(median(speed_sim{i},2,'omitnan'), median(speeds{i,j},2,'omitnan'), 'NMSE');
        nmse_omega(i,j) = goodnessOfFit(median(abs(omega_sim{i}(1:overlap,:)),2,'omitnan'), median(abs(omegas{i,j}(1:overlap,:)),2,'omitnan'), 'NMSE');
        nmse_total(i,j) = mean([nmse_speed(i,j), nmse_omega(i,j)]);
        
        mape_speed(i,j) = mape(median(speed_sim{i},2,'omitnan'), median(speeds{i,j},2,'omitnan'));
        mape_omega(i,j) = mape(median(abs(omega_sim{i}(1:overlap,:)),2,'omitnan'), median(abs(omegas{i,j}(1:overlap,:)),2,'omitnan'));
        mape_total(i,j) = mean([mape_speed(i,j), mape_omega(i,j)]);
        
        wmape_speed(i,j) = mape(median(speed_sim{i},2,'omitnan'), median(speeds{i,j},2,'omitnan'),'wMAPE');
        wmape_omega(i,j) = mape(median(abs(omega_sim{i}(1:overlap,:)),2,'omitnan'), median(abs(omegas{i,j}(1:overlap,:)),2,'omitnan'),'wMAPE');
        wmape_total(i,j) = mean([wmape_speed(i,j), wmape_omega(i,j)]);
    end
end

%% PRINT RESULTS
metrics_of_interest = {wmape_speed, wmape_total, wmape_omega}; metrics_tags = ["wmape_v", "wmape_{tot}", "wmape\omega"]; metrics_color = ['b','k','r'];

% % Single experiment plots
% for i = 1:size(experiments_names,1)  % for each experiment
if size(experiments_names,1) == 1
    fileID = fopen(fullfile(simulations_folder, sim_names(i), 'multi_exp_comparison.txt'),'wt');
    fprintf(fileID,'DOME multi experiment comparison\n\n');
    fprintf(fileID,'Date: %s\n',datestr(now, 'dd/mm/yy'));
    fprintf(fileID,'Time: %s\n\n',datestr(now, 'HH:MM'));
    fprintf(fileID,'Experiment\t\tNMSE speed\tNMSE omega\tNMSE tot\n');
    for j=1:size(experiments_names,2)
        fprintf(fileID,'%s\t',experiments_names(i,j));
        fprintf(fileID,'%.2f\t\t%.2f\t\t%.2f\n',nmse_speed(i,j),nmse_omega(i,j),nmse_total(i,j));
    end
    fclose(fileID);
    
    figure %time plot
    subplot(2,1,1)
    hold on
    ylim([0,100])
    highlightInputs(timeInstants, u{i}, 'r', 0.25)
    for j=2:size(experiments_names,2)
        plot(timeInstants, median(speeds{i,j},2,'omitnan'),'b', color=[0.5,0.5,1]);
    end
    l1=plot(timeInstants, median(speeds{i,1},2,'omitnan'),'b',LineWidth=2);
    l2=plot(timeInstants, median(speed_sim{i},2,'omitnan'),'k',LineWidth=2);
    xlim([0,max(timeInstants)])
    xlabel('$t$ [s]','Interpreter','Latex','FontSize',16)
    ylabel('$v$ [px/s]','Interpreter','Latex','FontSize',16)
    legend([l1,l2],'REAL','SIMULATED')
    box on
    subplot(2,1,2)
    hold on
    ylim([0,1.5])
    highlightInputs(timeInstants, u{i}, 'r', 0.25)
    for j=2:size(experiments_names,2)
        plot(timeInstants, median(abs(omegas{i,j}),2,'omitnan'),'b', color=[0.5,0.5,1]);
    end
    l1=plot(timeInstants, median(abs(omegas{i,1}),2,'omitnan'),'b',LineWidth=2);
    l2=plot(timeInstants, median(abs(omega_sim{i}),2,'omitnan'),'k',LineWidth=2);
    xlabel('$t$ [s]','Interpreter','Latex','FontSize',16)
    ylabel('$|\omega|$ [rad/s]','Interpreter','Latex','FontSize',16)
    legend([l1,l2],'REAL','SIMULATED')
    box on
    saveas(gcf,fullfile(simulations_folder,sim_names(i), 'multi_exp_comparison_time_plot'))
    saveas(gcf,fullfile(simulations_folder,sim_names(i), 'multi_exp_comparison_time_plot'),'png')
    
    
    figure % metrics scatter single exp
    hold on
    for k=1:length(metrics_of_interest)
        x_pos = [1-(length(metrics_of_interest)-1)*0.1+(k-1)*0.2]'-linspace(-1,1,size(experiments_names,2)-1)*0.033;
        plots(k,:)=bar(mean(x_pos,2),metrics_of_interest{k}(:,1),0.15,metrics_color(k),'FaceAlpha',0.5);
        scatter(x_pos,metrics_of_interest{k}(:,2:end),100,metrics_color(k),'MarkerFaceColor','w','LineWidth',1);
    end
    xticks([1])
    xticklabels(tags(i))
    set(gca, 'TickLabelInterpreter', 'none');
    xlim([0.5,1.5])
    all_metrics = [metrics_of_interest{:}];
    ylim([0, max(all_metrics(i,:),[],'all')*1.1])
    legend(plots(:,1),metrics_tags)
    set(gca,'FontSize',14)
    box on
    set(gca,'XGrid','off','YGrid','on')
    saveas(gcf,fullfile(fullfile(simulations_folder,sim_names(i)), 'multi_exp_comparison_NMSE'))
    saveas(gcf,fullfile(fullfile(simulations_folder,sim_names(i)), 'multi_exp_comparison_NMSE'),'png')
    
end

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
        legend([l1,l2],'REAL','SIMULATED')
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
        xlabel('$t$ [s]','Interpreter','Latex','FontSize',16)
        %ylabel('$|\omega|$ [rad/s]','Interpreter','Latex','FontSize',16)
        legend([l1,l2],'REAL','SIMULATED')
        box on
        
%         subplot(3,size(experiments_names,1),i+2*size(experiments_names,1))
%         hold on
%         for k=1:length(metrics_of_interest)
%             plots(k,:)=scatter(1-(length(metrics_of_interest)-1)*0.1+(k-1)*0.2,metrics_of_interest{k}(i,:),50,metrics_color(k));
%             scatter(1-(length(metrics_of_interest)-1)*0.1+(k-1)*0.2,metrics_of_interest{k}(i,1),50,metrics_color(k),"filled");
%         end
%         xticks([1])
%         xticklabels(tags(i))
%         set(gca, 'TickLabelInterpreter', 'none');
%         xlim([0,1+1])
%         ylim([0, max([metrics_of_interest{:}],[],'all')*1.1])
%         legend(plots(:,1),metrics_tags)
%         set(gca,'FontSize',14)
%         set(gca,'XGrid','off','YGrid','on')
%         box on
    end
    subplot(3,size(experiments_names,1),1); ylabel('$v$ [$\mu$m/s]','Interpreter','Latex','FontSize',14)
    subplot(3,size(experiments_names,1),1+size(experiments_names,1)); ylabel('$|\omega|$ [rad/s]','Interpreter','Latex','FontSize',14);
    subplot(3,size(experiments_names,1),[1+2*size(experiments_names,1),i+2*size(experiments_names,1)])
    hold on
    for k=1:length(metrics_of_interest)
        x_pos = [[1:length(tags)]-(length(metrics_of_interest)-1)*0.1+(k-1)*0.2]'-linspace(-1,1,size(experiments_names,2)-1)*0.025;
        plots(k,:)=bar(mean(x_pos,2),metrics_of_interest{k}(:,1),0.15,metrics_color(k),'FaceAlpha',0.5);
        scatter(x_pos,metrics_of_interest{k}(:,2:end),100,metrics_color(k),'MarkerFaceColor','w','LineWidth',1.25);
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
    saveas(gcf,fullfile(output_folder, 'multi_exp_comparison_overview'))
    saveas(gcf,fullfile(output_folder, 'multi_exp_comparison_overview'),'png')
    
end
