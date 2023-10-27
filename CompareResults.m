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

%% Add subfolders to the Matlab path
current_folder = fileparts(which('defaultParam'));
addpath(genpath(current_folder));

number_of_exp = length(experiment_paths);
experiments={};

%% Load data
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

figure % BOX PLOT - SPEED and ANGULAR VELOCITY - MEAN OVER AGENTS
subplot(2,1,1)
set(gca,'FontSize',14)
for i=1:number_of_exp; data_to_plot{i} = mean(experiments{i}.speed,'omitnan'); end
myboxplot(data_to_plot,true, 3)
ylabel('speed [px/s]')
xticklabels({'REAL','SIMULATED'})
title('Average over agents')
subplot(2,1,2)
set(gca,'FontSize',14)
for i=1:number_of_exp; data_to_plot{i} = mean(abs(experiments{i}.omega),'omitnan'); end
myboxplot(data_to_plot,true, 3)
ylabel('ang. vel. [rad/s]')
xticklabels({'REAL','SIMULATED'})
set(gcf,'position',[300,300,300,420])
if outputDir
    saveas(gcf,fullfile(outputDir, 'comparison_boxplot_meanOnAgents'))
    saveas(gcf,fullfile(outputDir, 'comparison_boxplot_meanOnAgents'),'png')
end

figure % BOX PLOT - SPEED and ANGULAR VELOCITY - MEAN OVER TIME
subplot(2,1,1)
set(gca,'FontSize',14)
for i=1:number_of_exp; data_to_plot{i} = mean(experiments{i}.speed,2,'omitnan'); end
myboxplot(data_to_plot,true, 3)
ylabel('speed [px/s]')
xticklabels({'REAL','SIMULATED'})
title('Average over time')
subplot(2,1,2)
set(gca,'FontSize',14)
for i=1:number_of_exp; data_to_plot{i} = mean(abs(experiments{i}.omega),2,'omitnan'); end
myboxplot(data_to_plot,true, 3)
ylabel('ang. vel. [rad/s]')
xticklabels({'REAL','SIMULATED'})
set(gcf,'position',[300,300,300,420])
if outputDir
    saveas(gcf,fullfile(outputDir, 'comparison_boxplot_meanOnTime'))
    saveas(gcf,fullfile(outputDir, 'comparison_boxplot_meanOnTime'),'png')
end

figure % BOX PLOT - SPEED and ANGULAR VELOCITY - ALL POINTS
subplot(2,1,1)
set(gca,'FontSize',14)
for i=1:number_of_exp; data_to_plot{i} = experiments{i}.speed(:); end
myboxplot(data_to_plot, true, 3)
ylabel('speed [px/s]')
xticklabels({'REAL','SIMULATED'})
title('All points')
subplot(2,1,2)
set(gca,'FontSize',14)
for i=1:number_of_exp; data_to_plot{i} = abs(experiments{i}.omega(:)); end
myboxplot(data_to_plot, true, 3)%, [0,0.4470,0.7410])
ylabel('ang. vel. [rad/s]')
xticklabels({'REAL','SIMULATED'})
set(gcf,'position',[300,300,300,420])
if outputDir
    saveas(gcf,fullfile(outputDir, 'comparison_boxplot_all'))
    saveas(gcf,fullfile(outputDir, 'comparison_boxplot_all'),'png')
end

figure % SCATTER PLOT - SPEED and ANGULAR VELOCITY - MEAN OVER TIME
x_to_plot=[];
y_to_plot=[];
grouping=[];
for i=1:number_of_exp; x_to_plot = [x_to_plot,  mean(experiments{i}.speed,2,'omitnan')']; end
for i=1:number_of_exp; y_to_plot = [y_to_plot,  mean(abs(experiments{i}.omega),2,'omitnan')']; end
for i=1:number_of_exp; grouping = [grouping,  i*ones(1,length(mean(experiments{i}.speed,2,'omitnan')))]; end
s=scatterhist(x_to_plot, y_to_plot, 'Location','NorthEast','Direction','out','Group', grouping ,'Kernel','on');%,'Marker','.', 'MarkerSize', 15);
xlabel(s,'mean speed [px/s]')
ylabel(s,'mean ang. vel. [rad/s]')
s(1).YAxisLocation = 'left';
s(1).XAxisLocation = 'bottom';
%set(get(gca,'children'),'filled',true)
s(2).Position = [0.1    0.82   0.7    0.125];
s(3).Position = [0.82   0.1    0.125    0.7];
s(1).Position(3) = 0.7;
s(1).Position(4) = 0.7;
legend({'REAL','SIMULATED'})
if outputDir
saveas(gcf,fullfile(outputDir, 'comparison_scatter_meanOnTime'))
saveas(gcf,fullfile(outputDir, 'comparison_scatter_meanOnTime'),'png')
end



