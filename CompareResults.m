%% COMPARE EXPERIMENTS

close all
clear

experiment_paths = {'/Volumes/DOMEPEN/Experiments/2023_06_15_Euglena_1/tracking_2023_10_12';
                '/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/2023_10_25_IndependentSDEs_4';
                '/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/2023_10_25_IndependentSDEs_3';
                '/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/2023_10_25_IndependentSDEs_2';
                '/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/2023_10_25_IndependentSDEs_1'
                };

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
myboxplot(data_to_plot)
ylabel('speed [px/s]')
xticklabels({'REAL','SIMULATED'})
subplot(2,1,2)
set(gca,'FontSize',14)
for i=1:number_of_exp; data_to_plot{i} = mean(abs(experiments{i}.omega),'omitnan'); end
myboxplot(data_to_plot)
ylabel('ang. vel. [rad/s]')
xticklabels({'REAL','SIMULATED'})

