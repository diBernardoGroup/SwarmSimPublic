clear
close all

simulations_folder = '/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations';
experiments_folder = "/Volumes/DOMEPEN/Experiments";

% switch 10s all
% tag = 'Euglena_switch_10';
% sim_name = '2024_04_17_switch_10_2';
% experiments_names = ["comparisons/Euglena_switch_10/combo5", "2023_06_15_Euglena_7", "2023_06_26_Euglena_23","2023_06_26_Euglena_24", "2023_07_10_Euglena_15", "2023_07_10_Euglena_16"];

% switch 10s selected
tag = 'switch_10';
sim_name = '2024_04_17_switch_10_2';
experiments_names = ["comparisons/Euglena_switch_10/combo3","2023_06_15_Euglena_7", "2023_06_26_Euglena_24", "2023_07_10_Euglena_15"];

% switch 5s
tag = 'Euglena_switch_5';
sim_name = '2024_04_17_switch_5_1';
experiments_names = ["comparisons/Euglena_switch_5/combo", "2023_06_15_Euglena_8", "2023_06_26_Euglena_25", "2023_07_10_Euglena_18"];

% switch 1s
tag = 'Euglena_switch_1';
sim_name = '2024_04_17_switch_1_2';
experiments_names = ["comparisons/Euglena_switch_1/combo", "2023_06_15_Euglena_11", "2023_06_26_Euglena_28", "2023_07_10_Euglena_19"];

% off
tag = 'Euglena_OFF';
sim_name = '2024_04_17_OFF_2';
experiments_names = ["comparisons/Euglena_OFF/combo", "2023_06_15_Euglena_1", "2023_06_26_Euglena_13", "2023_07_10_Euglena_6"];

% OFF-ON-OFF 75
tag = 'Euglena_75_ON';
sim_name = '2024_04_17_75_ON_1';
experiments_names = ["comparisons/Euglena_75_ON/combo", "2023_06_15_Euglena_2", "2023_06_26_Euglena_16", "2023_07_10_Euglena_8"];

% OFF-ON-OFF 150
tag = 'Euglena_150_ON';
sim_name = '2024_04_17_150_ON_3';
experiments_names = ["comparisons/Euglena_150_ON/combo", "2023_06_15_Euglena_3", "2023_06_26_Euglena_18", "2023_07_10_Euglena_10"];

% OFF-ON-OFF 255
tag = 'Euglena_255_ON';
sim_name = '2024_04_17_255_ON_3';
experiments_names = ["comparisons/Euglena_255_ON/combo", "2023_06_15_Euglena_4", "2023_06_26_Euglena_20", "2023_07_10_Euglena_12"];

% Ramp
tag = 'Euglena_ramp';
sim_name = '2024_04_17_ramp_2';
experiments_names = ["comparisons/Euglena_ramp/combo", "2023_06_15_Euglena_6", "2023_06_26_Euglena_22", "2023_07_10_Euglena_14"];

deltaT = 0.5;
timeInstants = [0:deltaT:180];

sim_folder = fullfile(simulations_folder,sim_name);
output_folder = sim_folder;

%% LOAD DATA
% for i=1:length(data_folders)   % for experiment
%     [NMSE_speed,NMSE_omega,NMSE_total] = compareResults({data_folders(i),sim_folder}, char(data_folders(1)));
% end

for i=1:length(experiments_names)   % for experiment
    if i==1
    data_folders(i) =  fullfile(experiments_folder,experiments_names(i));
    else
    subFolfder = getLastTracking(fullfile(experiments_folder,experiments_names(i)));
    data_folders(i) =  fullfile(experiments_folder,experiments_names(i),subFolfder);
    end
    speeds{i}  = load(fullfile(data_folders(i),'speeds_smooth.txt'));
    omegas{i}  = load(fullfile(data_folders(i),'ang_vel_smooth.txt'));
    
    if i==1
        inputs = load(fullfile(data_folders(i),'inputs.txt'));
        u=inputs(:,1)/255;              %select blue channel and scale in [0,1]
    else
        assert( all(inputs==load(fullfile(data_folders(i),'inputs.txt')),'all') )
    end
end


sim_data = load(fullfile(sim_folder,'data.mat'));
speed_sim = sim_data.speed;
omega_sim = sim_data.omega(1:end-1,:);

for i=1:length(data_folders)   % for experiment
    NMSE_speed(i) = goodnessOfFit(median(speed_sim,2,'omitnan'), median(speeds{i},2,'omitnan'), 'NMSE');
    NMSE_omega(i) = goodnessOfFit(median(abs(omega_sim),2,'omitnan'), median(abs(omegas{i}),2,'omitnan'), 'NMSE');
    NMSE_total(i) = mean([NMSE_speed(i), NMSE_omega(i)]);
    disp(['NMSE between median for speed:',num2str(NMSE_speed(i),'%.2f'),' and omega: ',num2str(NMSE_omega(i),'%.2f'),' total: ', num2str(NMSE_total(i),'%.2f')])
end


%% PRINT RESULTS

metrics_of_interest = {NMSE_speed, NMSE_omega, NMSE_total};
combine_metrics = true;
metrics_color = ['b','r','k'];
metrics_tags = ["NMSE_v", "NMSE_\omega", "NMSE_{tot}"];


fileID = fopen(fullfile(output_folder, 'multi_exp_comparison.txt'),'wt');
fprintf(fileID,'DOME multi experiment comparison\n\n');
fprintf(fileID,'Date: %s\n',datestr(now, 'dd/mm/yy'));
fprintf(fileID,'Time: %s\n\n',datestr(now, 'HH:MM'));
fprintf(fileID,'Experiment\t\t\t\t\tNMSE speed\tNMSE omega\tNMSE tot\n');
for i=1:length(experiments_names)
    fprintf(fileID,'%s\t',experiments_names(i));
    fprintf(fileID,'%.2f\t\t%.2f\t\t%.2f\n',NMSE_speed(i),NMSE_omega(i),NMSE_total(i));
end
fclose(fileID);

figure %time plot
subplot(2,1,1)
hold on
ylim([0,100])
if isvarname('u')
    highlightInputs(timeInstants, u, 'r', 0.25)
end
for i=2:length(experiments_names)
    plot(timeInstants, median(speeds{i},2,'omitnan'),'b', color=[0.5,0.5,1]);
end
l1=plot(timeInstants, median(speeds{1},2,'omitnan'),'b',LineWidth=2);
l2=plot(timeInstants, median(speed_sim,2,'omitnan'),'k',LineWidth=2);
xlim([0,max(timeInstants)])
xlabel('$t$ [s]','Interpreter','Latex','FontSize',16)
ylabel('$v$ [px/s]','Interpreter','Latex','FontSize',16)
legend([l1,l2],'REAL','SIMULATED')
box on
subplot(2,1,2)
hold on
ylim([0,1.5])
if isvarname('u')
    highlightInputs(timeInstants, u, 'r', 0.25)
end
for i=2:length(data_folders)
    plot(timeInstants(1:end-1), median(abs(omegas{i}),2,'omitnan'),'b', color=[0.5,0.5,1]);
end
l1=plot(timeInstants(1:end-1), median(abs(omegas{1}),2,'omitnan'),'b',LineWidth=2);
l2=plot(timeInstants(1:end-1), median(abs(omega_sim),2,'omitnan'),'k',LineWidth=2);
xlabel('$t$ [s]','Interpreter','Latex','FontSize',16)
ylabel('$|\omega|$ [rad/s]','Interpreter','Latex','FontSize',16)
legend([l1,l2],'REAL','SIMULATED')
box on
saveas(gcf,fullfile(output_folder, 'multi_exp_comparison_time_plot'))
saveas(gcf,fullfile(output_folder, 'multi_exp_comparison_time_plot'),'png')


figure % NMSE
hold on
for i=1:length(metrics_of_interest)
    plots(i,:)=scatter(1-(length(metrics_of_interest)-1)*0.1+(i-1)*0.2,metrics_of_interest{i},50,metrics_color(i));
    scatter(1-(length(metrics_of_interest)-1)*0.1+(i-1)*0.2,metrics_of_interest{i}(1),50,metrics_color(i),"filled");
end
xticks([1])
xticklabels(tag)
set(gca, 'TickLabelInterpreter', 'none');
xlim([0,1+1])
ylim([0, max([metrics_of_interest{:}],[],'all')*1.1])
legend(plots(:,1),metrics_tags)
set(gca,'FontSize',14)
box on
saveas(gcf,fullfile(output_folder, 'multi_exp_comparison_NMSE'))
saveas(gcf,fullfile(output_folder, 'multi_exp_comparison_NMSE'),'png')


