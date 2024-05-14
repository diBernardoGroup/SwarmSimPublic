clear
close all


data_folder = '/Volumes/DOMEPEN/Experiments/2023_06_15_Euglena_1/tracking_2023_10_12'; % off
%data_folder = '/Volumes/DOMEPEN/Experiments/2023_06_15_Euglena_7/tracking_2023_10_16'; % switch10s
data_folder = '/Volumes/DOMEPEN/Experiments/comparisons/Euglena_switch_10/combo5'; % switch10s combo
% data_folder = '/Volumes/DOMEPEN/Experiments/2023_06_26_Euglena_19/tracking_2023_10_16'; % on255

deltaT = 0.5;
dT = 0.01;

current_folder = fileparts(which('AnalyseDOMEexp'));
addpath(genpath(current_folder));

%% Load data
speed  = load(fullfile(data_folder,'speeds_smooth.txt'));
omega  = load(fullfile(data_folder,'ang_vel_smooth.txt'));
inputs = load(fullfile(data_folder,'inputs.txt'));

N = size(speed,2);
timeInstants = [0:size(speed,1)-1] * deltaT;
agents = [0:N-1]';
u=inputs(:,1)/255;
u_dot_BE = [0;diff(u)]/deltaT;
u_dot_grad = gradient(u)/deltaT;
u_dot = u_dot_BE;
u_dot = max(u_dot,0);
u_matrix = [u, u_dot];

%% PLOTS
outputDir = fullfile(data_folder,'plots');
if ~exist(outputDir,'dir'); mkdir(outputDir); end

figure % SCATTER PLOT - SPEED and ANGULAR VELOCITY - MEAN OVER TIME
x_to_plot=[];
y_to_plot=[];
grouping=[];
s=scatterhist(mean(speed,1,'omitnan'), mean(abs(omega),1,'omitnan'), 'Location','NorthEast','Direction','out' ,'Color','bk','Kernel','on');
xlabel(s,'mean $v$ [$\mu$m/s]','Interpreter','Latex','FontSize',16)
ylabel(s,'mean $|\omega|$ [rad/s]','Interpreter','Latex','FontSize',16)
s(1).YAxisLocation = 'left';
s(1).XAxisLocation = 'bottom';
%set(get(gca,'children'),'filled',true)
s(2).Position = [0.1    0.82   0.7    0.125];
s(3).Position = [0.82   0.1    0.125    0.7];
s(1).Position(3) = 0.7;
s(1).Position(4) = 0.7;
ylim([0,3])
xlim([0,120])
if outputDir
    saveas(gcf,fullfile(outputDir, 'scatter_meanOnTime'))
    saveas(gcf,fullfile(outputDir, 'scatter_meanOnTime'),'png')
end

figure % TIME PLOT - SPEED and ANGULAR VELOCITY
subplot(2,1,1)
xlim([0,max(timeInstants)])
ylim([0,120])
if isvarname('u')
    highlightInputs(timeInstants, u, 'r', 0.25)
end
l1=plotWithShade(timeInstants, median(speed,2,'omitnan'), quantile(speed, 0.1, 2), quantile(speed, 0.9, 2), 'b', 0.3);
xlabel('$t$ [s]','Interpreter','Latex','FontSize',16)
ylabel('$v$ [$\mu$m/s]','Interpreter','Latex','FontSize',16)
rng=ylim;
box on
subplot(2,1,2)
xlim([0,max(timeInstants)])
ylim([0,2])
if isvarname('u')
    highlightInputs(timeInstants, u, 'r', 0.25)
end
l1=plotWithShade(timeInstants, median(abs(omega),2,'omitnan'), quantile(abs(omega), 0.1, 2), quantile(abs(omega), 0.9, 2), 'b', 0.3);
xlabel('$t$ [s]','Interpreter','Latex','FontSize',16)
ylabel('$|\omega|$ [rad/s]','Interpreter','Latex','FontSize',16)
rng=ylim;
box on
if outputDir
    saveas(gcf,fullfile(outputDir, 'time_plot'))
    saveas(gcf,fullfile(outputDir, 'time_plot'),'png')
end

% figure % inputs
% subplot(3,1,1)
% hold on
% plot(timeInstants,u,'--')
% plot(t_sim,u_sim(:,1),'--')
% legend({'experiment','simulation'})
% subplot(3,1,2)
% hold on
% plot(timeInstants,u_dot,'--')
% plot(t_sim,u_sim(:,2),'--')
% legend({'experiment','simulation'})
% subplot(3,1,3)
% hold on
% plot(timeInstants,cumtrapz(u_dot)*deltaT,'--')
% plot(t_sim,cumtrapz(u_sim(:,2))*dT,'--')
% legend({'experiment','simulation'})
