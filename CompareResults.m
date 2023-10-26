%% COMPARE EXPERIMENTS

close all
clear

experiment_paths = {'/Volumes/DOMEPEN/Experiments/2023_06_15_Euglena_1/tracking_2023_10_12';
                '/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/2023_10_26_IndependentSDEs_9';
                %'/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/2023_10_26_IndependentSDEs_8';
                %'/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/2023_10_26_IndependentSDEs_7';
                %'/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/2023_10_26_IndependentSDEs_6';
                %'/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/2023_10_26_IndependentSDEs_5';
                %'/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/2023_10_26_IndependentSDEs_4';
                %'/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/2023_10_26_IndependentSDEs_3';
                %'/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/2023_10_26_IndependentSDEs_1';
                %'/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/2023_10_25_IndependentSDEs_4';
                %'/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/2023_10_25_IndependentSDEs_3';
                %'/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/2023_10_25_IndependentSDEs_2';
                %'/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/2023_10_25_IndependentSDEs_1'
                };
            
outputDir='/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/comparison';

number_of_exp = length(experiment_paths);
experiments={};

% load data
for i=1:number_of_exp
    exp = experiment_paths{i};
    data_path = fullfile(exp,'data.mat');
    if isfile(data_path)
        d = load(data_path);
        omega = d.omega';
        speed = d.speed';
    else
        speed   = load(fullfile(exp,'speeds_smooth.txt'));
        speed   = speed(1:end-1,:)';
        omega = load(fullfile(exp,'ang_vel_smooth.txt'))';
    end
    
    data = table(speed, omega);
    data(all(isnan(data.speed),2),:) = [];
    experiments{i}=data;
end

figure
subplot(2,1,1)
set(gca,'FontSize',14)
for i=1:number_of_exp; data_to_plot{i} = mean(experiments{i}.speed,'omitnan'); end
myboxplot(data_to_plot,true, 3, [0,0.4470,0.7410])
ylabel('speed [px/s]')
xticklabels({'REAL','SIMULATED'})
title('Average over agents')
subplot(2,1,2)
set(gca,'FontSize',14)
for i=1:number_of_exp; data_to_plot{i} = mean(abs(experiments{i}.omega),'omitnan'); end
myboxplot(data_to_plot,true, 3, [0,0.4470,0.7410])
ylabel('ang. vel. [rad/s]')
xticklabels({'REAL','SIMULATED'})
set(gcf,'position',[300,300,300,420])
if outputDir
    saveas(gcf,fullfile(outputDir, 'comparison_boxplot_meanOnAgents'))
    saveas(gcf,fullfile(outputDir, 'comparison_boxplot_meanOnAgents'),'png')
end

figure
subplot(2,1,1)
set(gca,'FontSize',14)
for i=1:number_of_exp; data_to_plot{i} = mean(experiments{i}.speed,2,'omitnan'); end
myboxplot(data_to_plot,true, 3, [0,0.4470,0.7410])
ylabel('speed [px/s]')
xticklabels({'REAL','SIMULATED'})
title('Average over time')
subplot(2,1,2)
set(gca,'FontSize',14)
for i=1:number_of_exp; data_to_plot{i} = mean(abs(experiments{i}.omega),2,'omitnan'); end
myboxplot(data_to_plot,true, 3, [0,0.4470,0.7410])
ylabel('ang. vel. [rad/s]')
xticklabels({'REAL','SIMULATED'})
set(gcf,'position',[300,300,300,420])
if outputDir
    saveas(gcf,fullfile(outputDir, 'comparison_boxplot_meanOnTime'))
    saveas(gcf,fullfile(outputDir, 'comparison_boxplot_meanOnTime'),'png')
end

figure
subplot(2,1,1)
set(gca,'FontSize',14)
for i=1:number_of_exp; data_to_plot{i} = experiments{i}.speed(:); end
myboxplot(data_to_plot, true, 3,[0,0.4470,0.7410])
ylabel('speed [px/s]')
xticklabels({'REAL','SIMULATED'})
title('All points')
subplot(2,1,2)
set(gca,'FontSize',14)
for i=1:number_of_exp; data_to_plot{i} = abs(experiments{i}.omega(:)); end
myboxplot(data_to_plot, true, 3, [0,0.4470,0.7410])
ylabel('ang. vel. [rad/s]')
xticklabels({'REAL','SIMULATED'})
set(gcf,'position',[300,300,300,420])
if outputDir
    saveas(gcf,fullfile(outputDir, 'comparison_boxplot_all'))
    saveas(gcf,fullfile(outputDir, 'comparison_boxplot_all'),'png')
end

