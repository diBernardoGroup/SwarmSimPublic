%% COMPARE EXPERIMENTS

function [mse_speed,mse_omega,nmse_speed,nmse_omega,nmse_total] = compareResults(experiment_paths, outputDir)

% experiment_paths = {'/Volumes/DOMEPEN/Experiments/2023_06_15_Euglena_1/tracking_2023_10_12';
%                 '/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/2023_10_26_IndependentSDEs_9';
%                 '/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/2023_10_26_IndependentSDEs_8';
%                 '/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/2023_10_26_IndependentSDEs_7';
%                 %'/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/2023_10_26_IndependentSDEs_6';
%                 %'/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/2023_10_26_IndependentSDEs_5';
%                 %'/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/2023_10_26_IndependentSDEs_4';
%                 %'/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/2023_10_26_IndependentSDEs_3';
%                 %'/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/2023_10_26_IndependentSDEs_1';
%                 %'/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/2023_10_25_IndependentSDEs_4';
%                 %'/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/2023_10_25_IndependentSDEs_3';
%                 %'/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/2023_10_25_IndependentSDEs_2';
%                 %'/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/2023_10_25_IndependentSDEs_1'
%                 };
            
% outputDir='/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/comparison';

arguments
    experiment_paths
    outputDir = ''
end

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
        timeInstants = d.timeInstants;
        u = d.u;
    else
        speed = load(fullfile(exp,'speeds_smooth.txt'));
        speed = speed(1:end-1,:)';
        omega = load(fullfile(exp,'ang_vel_smooth.txt'))';
    end
    
    data = table(speed, omega);
    data(all(isnan(data.speed),2),:) = [];
    experiments{i}=data;
end

if number_of_exp>2
    mse_speed = nan(number_of_exp);
    mse_omega = nan(number_of_exp);
    for i=1:number_of_exp
        for j=1:number_of_exp
            overlap = min(size(experiments{i}.speed,2),size(experiments{j}.speed,2));
            mse_speed(i,j) = mean((mean(experiments{1}.speed(:,1:overlap),1,'omitnan')-mean(experiments{2}.speed(:,1:overlap),1,'omitnan')).^2);
            mse_omega(i,j) = mean((mean(experiments{1}.omega(:,1:overlap),1,'omitnan')-mean(experiments{2}.omega(:,1:overlap),1,'omitnan')).^2); 
        end
    end
else
    overlap = min(size(experiments{1}.speed,2),size(experiments{2}.speed,2));
    mse_speed = mean((mean(experiments{1}.speed(:,1:overlap),1,'omitnan')-mean(experiments{2}.speed(:,1:overlap),1,'omitnan')).^2);
    mse_omega = mean((mean(experiments{1}.omega(:,1:overlap),1,'omitnan')-mean(experiments{2}.omega(:,1:overlap),1,'omitnan')).^2);
    nmse_speed = goodnessOfFit(mean(experiments{2}.speed(:,1:overlap),1,'omitnan')', mean(experiments{1}.speed(:,1:overlap),1,'omitnan')', 'NMSE');
    nmse_omega = goodnessOfFit(mean(experiments{2}.omega(:,1:overlap),1,'omitnan')', mean(experiments{1}.omega(:,1:overlap),1,'omitnan')', 'NMSE');
    nmse_total = mean([nmse_speed, nmse_omega]);
    %disp(['MSE between mean for speed:',num2str(mse_speed),' and omega: ',num2str(mse_omega)])
    disp(['NMSE between mean for speed:',num2str(nmse_speed),' and omega: ',num2str(nmse_omega),' total: ', num2str(nmse_total)])
    
    nmse_speed_med = goodnessOfFit(median(experiments{2}.speed(:,1:overlap),1,'omitnan')', median(experiments{1}.speed(:,1:overlap),1,'omitnan')', 'NMSE');
    nmse_omega_med = goodnessOfFit(median(experiments{2}.omega(:,1:overlap),1,'omitnan')', median(experiments{1}.omega(:,1:overlap),1,'omitnan')', 'NMSE');
    nmse_total_med = mean([nmse_speed_med, nmse_omega_med]);
    disp(['NMSE between median for speed:',num2str(nmse_speed_med),' and omega: ',num2str(nmse_omega_med),' total: ', num2str(nmse_total_med)])
end


figure % BOX PLOT - SPEED and ANGULAR VELOCITY - MEAN OVER AGENTS
subplot(2,1,1)
set(gca,'FontSize',12)
for i=1:number_of_exp; data_to_plot{i} = mean(experiments{i}.speed,'omitnan'); end
myboxplot(data_to_plot,true, 3, {'b','k'})
ylabel('$v$ [px/s]','Interpreter','Latex','FontSize',16)
xticklabels({'REAL','SIMULATED'})
title('Average over agents')
subplot(2,1,2)
set(gca,'FontSize',12)
for i=1:number_of_exp; data_to_plot{i} = mean(abs(experiments{i}.omega),'omitnan'); end
myboxplot(data_to_plot,true, 3, {'b','k'})
ylabel('$\omega$ [rad/s]','Interpreter','Latex','FontSize',16)
xticklabels({'REAL','SIMULATED'})
set(gcf,'position',[300,300,300,420])
if outputDir
    saveas(gcf,fullfile(outputDir, 'comparison_boxplot_meanOnAgents'))
    saveas(gcf,fullfile(outputDir, 'comparison_boxplot_meanOnAgents'),'png')
end

figure % BOX PLOT - SPEED and ANGULAR VELOCITY - MEAN OVER TIME
subplot(2,1,1)
set(gca,'FontSize',12)
for i=1:number_of_exp; data_to_plot{i} = mean(experiments{i}.speed,2,'omitnan'); end
myboxplot(data_to_plot,true, 3, {'b','k'})
ylabel('$v$ [px/s]','Interpreter','Latex','FontSize',16)
xticklabels({'REAL','SIMULATED'})
title('Average over time')
subplot(2,1,2)
set(gca,'FontSize',12)
for i=1:number_of_exp; data_to_plot{i} = mean(abs(experiments{i}.omega),2,'omitnan'); end
myboxplot(data_to_plot,true, 3, {'b','k'})
ylabel('$\omega$ [rad/s]','Interpreter','Latex','FontSize',16)
xticklabels({'REAL','SIMULATED'})
set(gcf,'position',[300,300,300,420])
if outputDir
    saveas(gcf,fullfile(outputDir, 'comparison_boxplot_meanOnTime'))
    saveas(gcf,fullfile(outputDir, 'comparison_boxplot_meanOnTime'),'png')
end

figure % BOX PLOT - SPEED and ANGULAR VELOCITY - ALL POINTS
subplot(2,1,1)
set(gca,'FontSize',12)
for i=1:number_of_exp; data_to_plot{i} = experiments{i}.speed(:); end
myboxplot(data_to_plot, true, 3, {'b','k'})
ylabel('$v$ [px/s]','Interpreter','Latex','FontSize',16)
xticklabels({'REAL','SIMULATED'})
title('All points')
subplot(2,1,2)
set(gca,'FontSize',12)
for i=1:number_of_exp; data_to_plot{i} = abs(experiments{i}.omega(:)); end
myboxplot(data_to_plot, true, 3, {'b','k'})%, [0,0.4470,0.7410])
ylabel('$\omega$ [rad/s]','Interpreter','Latex','FontSize',16)
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
s=scatterhist(x_to_plot, y_to_plot, 'Location','NorthEast','Direction','out','Group', grouping ,'Color','bk','Kernel','on');%,'Marker','.', 'MarkerSize', 15);
xlabel(s,'mean $v$ [px/s]','Interpreter','Latex','FontSize',16)
ylabel(s,'mean $\omega$ [rad/s]','Interpreter','Latex','FontSize',16)
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

figure % TIME PLOT - SPEED and ANGULAR VELOCITY
subplot(2,1,1)
l1=plotWithShade(timeInstants, median(experiments{1}.speed,1,'omitnan'), min(experiments{1}.speed, [], 1,'omitnan'), max(experiments{1}.speed, [], 1,'omitnan'), 'b', 0.3);
l2=plotWithShade(timeInstants, median(experiments{2}.speed,1,'omitnan'), min(experiments{2}.speed, [], 1,'omitnan'), max(experiments{2}.speed, [], 1,'omitnan'), 'k', 0.3);
xlim([0,max(timeInstants)])
if isvarname('u')
    highlightInputs(timeInstants, u, 'r', 0.25)
end
xlabel('$t$ [s]','Interpreter','Latex','FontSize',16)
ylabel('$v$ [px/s]','Interpreter','Latex','FontSize',16)
legend([l1,l2],'REAL','SIMULATED')
rng=ylim;
box on
% subplot(2,4,4)
% hold on
% h=histogram(experiments{1}.speed(:),'Orientation','horizontal');
% h=histogram(experiments{2}.speed(:),'Orientation','horizontal');
% ylim(rng);
% set(gca,'xtick',[])
subplot(2,1,2)
l1=plotWithShade(timeInstants, median(abs(experiments{1}.omega),1,'omitnan'), min(abs(experiments{1}.omega), [], 1,'omitnan'), max(abs(experiments{1}.omega), [], 1,'omitnan'), 'b', 0.3);
l2=plotWithShade(timeInstants, median(abs(experiments{2}.omega),1,'omitnan'), min(abs(experiments{2}.omega), [], 1,'omitnan'), max(abs(experiments{2}.omega), [], 1,'omitnan'), 'k', 0.3);
xlim([0,max(timeInstants)])
if isvarname('u')
    highlightInputs(timeInstants, u, 'r', 0.25)
end
xlabel('$t$ [s]','Interpreter','Latex','FontSize',16)
ylabel('$\omega$ [rad/s]','Interpreter','Latex','FontSize',16)
legend([l1,l2],'REAL','SIMULATED')
rng=ylim;
box on
% subplot(2,4,8)
% hold on
% h=histogram(abs(experiments{1}.omega(:)),'Orientation','horizontal', 'FaceColor', 'b', 'EdgeColor', 'none');
% h=histogram(abs(experiments{2}.omega(:)),'Orientation','horizontal', 'FaceColor', 'k', 'EdgeColor', 'none');
% ylim(rng);
% set(gca,'xtick',[])
if outputDir
    saveas(gcf,fullfile(outputDir, 'comparison_time_plot'))
    saveas(gcf,fullfile(outputDir, 'comparison_time_plot'),'png')
end

end


