clear
close all

simulations_folder = '/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations';
experiments_folder = "/Volumes/DOMEPEN/Experiments";

% switch 10s all
% tag = 'switch_10';
% sim_name = '2024_04_17_switch_10_2';
% experiments_names = ["comparisons/Euglena_switch_10/combo5", "2023_06_15_Euglena_7", "2023_06_26_Euglena_23","2023_06_26_Euglena_24", "2023_07_10_Euglena_15", "2023_07_10_Euglena_16"];

% % switch 10s selected
% tag = 'switch_10';
% sim_name = '2024_04_17_switch_10_2';
% experiments_names = ["comparisons/Euglena_switch_10/combo3","2023_06_15_Euglena_7", "2023_06_26_Euglena_24", "2023_07_10_Euglena_15"];
% 
% % switch 5s
% tag = 'switch_5';
% sim_name = '2024_04_17_switch_5_1';
% experiments_names = ["comparisons/Euglena_switch_5/combo", "2023_06_15_Euglena_8", "2023_06_26_Euglena_25", "2023_07_10_Euglena_18"];
% 
% % switch 1s
% tag = 'switch_1';
% sim_name = '2024_04_17_switch_1_2';
% experiments_names = ["comparisons/Euglena_switch_1/combo", "2023_06_15_Euglena_11", "2023_06_26_Euglena_28", "2023_07_10_Euglena_19"];
% 
% % off
% tag = 'OFF';
% sim_name = '2024_04_17_OFF_2';
% experiments_names = ["comparisons/Euglena_OFF/combo", "2023_06_15_Euglena_1", "2023_06_26_Euglena_13", "2023_07_10_Euglena_6"];
% 
% % OFF-ON-OFF 75
% tag = '75_ON';
% sim_name = '2024_04_17_75_ON_1';
% experiments_names = ["comparisons/Euglena_75_ON/combo", "2023_06_15_Euglena_2", "2023_06_26_Euglena_16", "2023_07_10_Euglena_8"];
% 
% % OFF-ON-OFF 150
% tag = '150_ON';
% sim_name = '2024_04_17_150_ON_3';
% experiments_names = ["comparisons/Euglena_150_ON/combo", "2023_06_15_Euglena_3", "2023_06_26_Euglena_18", "2023_07_10_Euglena_10"];
% 
% % OFF-ON-OFF 255
% tag = '255_ON';
% sim_name = '2024_04_17_255_ON_3';
% experiments_names = ["comparisons/Euglena_255_ON/combo", "2023_06_15_Euglena_4", "2023_06_26_Euglena_20", "2023_07_10_Euglena_12"];
% 
% % Ramp
% tag = 'ramp';
% sim_name = '2024_04_17_ramp_2';
%experiments_names = ["comparisons/Euglena_ramp/combo", "2023_06_15_Euglena_6", "2023_06_26_Euglena_22", "2023_07_10_Euglena_14"];

% all
tags = ["switch_10","switch_5"];
sim_names = ["2024_04_17_switch_10_2", "2024_04_17_switch_5_1"];
experiments_names = ["comparisons/Euglena_switch_10/combo3","2023_06_15_Euglena_7", "2023_06_26_Euglena_24", "2023_07_10_Euglena_15";
                     "comparisons/Euglena_switch_5/combo", "2023_06_15_Euglena_8", "2023_06_26_Euglena_25", "2023_07_10_Euglena_18"];

                 
output_folder = '/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/comparison';

deltaT = 0.5;
timeInstants = [0:deltaT:180];


%% LOAD DATA

for i = 1:size(experiments_names,1)  % for each experiment

% load simulation data
sim_folder = fullfile(simulations_folder,sim_names(i));
sim_data = load(fullfile(sim_folder,'data.mat'));
speed_sim{i} = sim_data.speed;
omega_sim{i} = sim_data.omega(1:end-1,:);

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
    
    if j==1
        inputs = load(fullfile(data_folder,'inputs.txt'));
        u{i}=inputs(:,1)/255;              %select blue channel and scale in [0,1]
    else
        assert( all(inputs==load(fullfile(data_folder,'inputs.txt')),'all') )
    end
end

% evaluate quality of fit
for j=1:size(experiments_names,2)   % for each replicate
    NMSE_speed(i,j) = goodnessOfFit(median(speed_sim{i},2,'omitnan'), median(speeds{i,j},2,'omitnan'), 'NMSE');
    NMSE_omega(i,j) = goodnessOfFit(median(abs(omega_sim{i}),2,'omitnan'), median(abs(omegas{i,j}),2,'omitnan'), 'NMSE');
    NMSE_total(i,j) = mean([NMSE_speed(i,j), NMSE_omega(i,j)]);
    %disp(['NMSE between median for speed:',num2str(NMSE_speed(i,j),'%.2f'),' and omega: ',num2str(NMSE_omega(i,j),'%.2f'),' total: ', num2str(NMSE_total(i,j),'%.2f')])
end
end

%% PRINT RESULTS

metrics_of_interest = {NMSE_speed, NMSE_omega, NMSE_total};
combine_metrics = true;
metrics_color = ['b','r','k'];
metrics_tags = ["NMSE_v", "NMSE_\omega", "NMSE_{tot}"];

% Single experiment plots
for i = 1:size(experiments_names,1)  % for each experiment

fileID = fopen(fullfile(simulations_folder, sim_names(i), 'multi_exp_comparison.txt'),'wt');
fprintf(fileID,'DOME multi experiment comparison\n\n');
fprintf(fileID,'Date: %s\n',datestr(now, 'dd/mm/yy'));
fprintf(fileID,'Time: %s\n\n',datestr(now, 'HH:MM'));
fprintf(fileID,'Experiment\t\tNMSE speed\tNMSE omega\tNMSE tot\n');
for j=1:size(experiments_names,2)
    fprintf(fileID,'%s\t',experiments_names(i,j));
    fprintf(fileID,'%.2f\t\t%.2f\t\t%.2f\n',NMSE_speed(i,j),NMSE_omega(i,j),NMSE_total(i,j));
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
    plot(timeInstants(1:end-1), median(abs(omegas{i,j}),2,'omitnan'),'b', color=[0.5,0.5,1]);
end
l1=plot(timeInstants(1:end-1), median(abs(omegas{i,1}),2,'omitnan'),'b',LineWidth=2);
l2=plot(timeInstants(1:end-1), median(abs(omega_sim{i}),2,'omitnan'),'k',LineWidth=2);
xlabel('$t$ [s]','Interpreter','Latex','FontSize',16)
ylabel('$|\omega|$ [rad/s]','Interpreter','Latex','FontSize',16)
legend([l1,l2],'REAL','SIMULATED')
box on
saveas(gcf,fullfile(simulations_folder,sim_names(i), 'multi_exp_comparison_time_plot'))
saveas(gcf,fullfile(simulations_folder,sim_names(i), 'multi_exp_comparison_time_plot'),'png')


figure % NMSE scatter single exp
hold on
for k=1:length(metrics_of_interest)
    plots(k,:)=scatter(1-(length(metrics_of_interest)-1)*0.1+(k-1)*0.2,metrics_of_interest{k}(i,:),50,metrics_color(k));
    scatter(1-(length(metrics_of_interest)-1)*0.1+(k-1)*0.2,metrics_of_interest{k}(i,1),50,metrics_color(k),"filled");
end
xticks([1])
xticklabels(tags(i))
set(gca, 'TickLabelInterpreter', 'none');
xlim([0,1+1])
all_metrics = [metrics_of_interest{:}];
ylim([0, max(all_metrics(i,:),[],'all')*1.1])
legend(plots(:,1),metrics_tags)
set(gca,'FontSize',14)
box on
set(gca,'XGrid','off','YGrid','on')
saveas(gcf,fullfile(fullfile(simulations_folder,sim_names(i)), 'multi_exp_comparison_NMSE'))
saveas(gcf,fullfile(fullfile(simulations_folder,sim_names(i)), 'multi_exp_comparison_NMSE'),'png')

end

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
xlabel('$t$ [s]','Interpreter','Latex','FontSize',16)
ylabel('$v$ [px/s]','Interpreter','Latex','FontSize',16)
legend([l1,l2],'REAL','SIMULATED')
box on
subplot(3,size(experiments_names,1),i+size(experiments_names,1))
hold on
ylim([0,1.5])
highlightInputs(timeInstants, u{i}, 'r', 0.25)
for j=2:size(experiments_names,2)
    plot(timeInstants(1:end-1), median(abs(omegas{i,j}),2,'omitnan'),'b', color=[0.5,0.5,1]);
end
l1=plot(timeInstants(1:end-1), median(abs(omegas{i,1}),2,'omitnan'),'b',LineWidth=2);
l2=plot(timeInstants(1:end-1), median(abs(omega_sim{i}),2,'omitnan'),'k',LineWidth=2);
xlabel('$t$ [s]','Interpreter','Latex','FontSize',16)
ylabel('$|\omega|$ [rad/s]','Interpreter','Latex','FontSize',16)
legend([l1,l2],'REAL','SIMULATED')
box on

subplot(3,size(experiments_names,1),i+2*size(experiments_names,1))
hold on
for k=1:length(metrics_of_interest)
    plots(k,:)=scatter(1-(length(metrics_of_interest)-1)*0.1+(k-1)*0.2,metrics_of_interest{k}(i,:),50,metrics_color(k));
    scatter(1-(length(metrics_of_interest)-1)*0.1+(k-1)*0.2,metrics_of_interest{k}(i,1),50,metrics_color(k),"filled");
end
xticks([1])
xticklabels(tags(i))
set(gca, 'TickLabelInterpreter', 'none');
xlim([0,1+1])
ylim([0, max([metrics_of_interest{:}],[],'all')*1.1])
legend(plots(:,1),metrics_tags)
set(gca,'FontSize',14)
set(gca,'XGrid','off','YGrid','on')
box on
end
saveas(gcf,fullfile(output_folder, 'multi_exp_comparison_overview'))
saveas(gcf,fullfile(output_folder, 'multi_exp_comparison_overview'),'png')


figure % NMSE scatter comparison
hold on
for k=1:length(metrics_of_interest)
    plots(k,:)=scatter([1:length(tags)]-(length(metrics_of_interest)-1)*0.1+(k-1)*0.2,metrics_of_interest{k},50,metrics_color(k));
    scatter([1:length(tags)]-(length(metrics_of_interest)-1)*0.1+(k-1)*0.2,metrics_of_interest{k}(:,1),50,metrics_color(k),"filled");
end
xticks([1:length(tags)])
xticklabels(tags)
set(gca, 'TickLabelInterpreter', 'none');
xlim([0,length(tags)+1])
ylim([0, max([metrics_of_interest{:}],[],'all')*1.1])
legend(plots(:,1),metrics_tags)
set(gca,'FontSize',14)
box on
set(gca,'XGrid','off','YGrid','on')
saveas(gcf,fullfile(output_folder, 'multi_exp_comparison_NMSE'))
saveas(gcf,fullfile(output_folder, 'multi_exp_comparison_NMSE'),'png')

