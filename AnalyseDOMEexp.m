% clear
% clc
% close all


current_folder = '/Volumes/DOMEPEN/Experiments/2023_06_15_Euglena_1/tracking_2023_10_12';

deltaT=0.5;

%% Load data
%identification=readtable(fullfile(current_folder,'identification.txt'));
load(fullfile(current_folder,'speeds_smooth.txt'))
load(fullfile(current_folder,'ang_vel_smooth.txt'))

N = size(speeds_smooth,2);
ids = [0:N-1]';

%% Identification
[mu_s, theta_s, sigma_s] = MLE_SDE_parameters(speeds_smooth, deltaT);
[mu_w, theta_w, sigma_w] = MLE_SDE_parameters(ang_vel_smooth, deltaT);
mu_s=mu_s'; theta_s=theta_s'; sigma_s=sigma_s'; mu_w=mu_w'; theta_w=theta_w'; sigma_w=sigma_w'; 
mu_s=round(mu_s,4); theta_s=round(theta_s,4); sigma_s=round(sigma_s,4); mu_w=round(mu_w,4); theta_w=round(theta_w,4); sigma_w=round(sigma_w,4) ; 

identification = table(ids, mu_s, theta_s, sigma_s, mu_w, theta_w, sigma_w);
identification(isnan(identification.mu_s),:) = [];
writetable(identification,fullfile(current_folder, 'identification.txt') ,'Delimiter',' ')
 

if exist('SDEparameters')
    figure
    subplot(2,1,1)
    set(gca,'FontSize',14)
    boxplot([[SDEparameters.mean_s-SDEparameters.std_s, SDEparameters.mean_s+SDEparameters.std_s]', [mean(speed)-std(speed);mean(speed)+std(speed)]])
    hold on
    xline(1.5)
    ylabel('speed [px/s]')
    xticks([1, size(speed,2)/2+1.5])
    xticklabels({'REAL','SIMULATED'})
    subplot(2,1,2)
    set(gca,'FontSize',14)
    boxplot([[SDEparameters.mean_w-SDEparameters.std_w, SDEparameters.mean_w+SDEparameters.std_w]', [mean(omega)-std(omega);mean(omega)+std(omega)]])
    hold on
    xline(1.5)
    ylabel('ang. vel. [rad/s]')
    xticks([1, size(speed,2)/2+1.5])
    xticklabels({'REAL','SIMULATED'})
end

