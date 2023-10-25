clear
close all


current_folder = '/Volumes/DOMEPEN/Experiments/2023_06_15_Euglena_1/tracking_2023_10_12';

deltaT=0.5;

%% Load data
%identification=readtable(fullfile(current_folder,'identification.txt'));
load(fullfile(current_folder,'speeds_smooth.txt'))
load(fullfile(current_folder,'ang_vel_smooth.txt'))

N = size(speeds_smooth,2);
agents = [0:N-1]';

%% Identification
[mu_s, theta_s, sigma_s] = LS_SDE_parameters(speeds_smooth, deltaT);
[mu_w, theta_w, sigma_w] = LS_SDE_parameters(ang_vel_smooth, deltaT);
mu_s=mu_s'; theta_s=theta_s'; sigma_s=sigma_s'; mu_w=mu_w'; theta_w=theta_w'; sigma_w=sigma_w'; 
mu_s=round(mu_s,4); theta_s=round(theta_s,4); sigma_s=round(sigma_s,4); mu_w=round(mu_w,4); theta_w=round(theta_w,4); sigma_w=round(sigma_w,4) ; 

mean_s  = round(mean(speeds_smooth,'omitnan')',4);
std_s   = round(std(speeds_smooth,'omitnan')',4);
mean_w  = round(mean(ang_vel_smooth,'omitnan')',4);
std_w   = round(std(ang_vel_smooth,'omitnan')',4);

identification = table(agents, mu_s, theta_s, sigma_s, mu_w, theta_w, sigma_w, mean_s, std_s, mean_w, std_w);
identification(isnan(identification.mu_s),:) = [];
writetable(identification,fullfile(current_folder, 'identification.txt') ,'Delimiter',' ')

identification