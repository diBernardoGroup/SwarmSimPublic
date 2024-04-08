clear
close all

sim_folder = '/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/2024_04_08_IndependentSDEsWithInput_8'; % folder simulation data
data_folders = ["/Volumes/DOMEPEN/Experiments/comparisons/Euglena_switch_10/combo",
    "/Volumes/DOMEPEN/Experiments/2023_06_15_Euglena_7/tracking_2023_10_16",
    "/Volumes/DOMEPEN/Experiments/2023_06_26_Euglena_23/tracking_2023_10_12",
    "/Volumes/DOMEPEN/Experiments/2023_06_26_Euglena_24/tracking_2024_04_08",
    "/Volumes/DOMEPEN/Experiments/2023_07_10_Euglena_15/tracking_2023_10_12",
    "/Volumes/DOMEPEN/Experiments/2023_07_10_Euglena_16/tracking_2024_04_08"]; % experiments data folders

dT = 0.01;
deltaT = 0.5;
timeInstants = [0:deltaT:180];

%% LOAD DATA
% for i=1:length(data_folders)   % for experiment
%     [NMSE_speed,NMSE_omega,NMSE_total] = compareResults({data_folders(i),sim_folder}, char(data_folders(1)));
% end

for i=1:length(data_folders)   % for experiment
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
    nmse_speed_med(i) = goodnessOfFit(median(speed_sim,2,'omitnan'), median(speeds{i},2,'omitnan'), 'NMSE');
    nmse_omega_med(i) = goodnessOfFit(median(abs(omega_sim),2,'omitnan'), median(abs(omegas{i}),2,'omitnan'), 'NMSE');
    nmse_total_med(i) = mean([nmse_speed_med(i), nmse_omega_med(i)]);
    disp(['NMSE between median for speed:',num2str(nmse_speed_med(i),'%.2f'),' and omega: ',num2str(nmse_omega_med(i),'%.2f'),' total: ', num2str(nmse_total_med(i),'%.2f')])
end




%% PRINT RESULTS

figure %time plot
subplot(2,1,1)
hold on
for i=2:length(data_folders)
    plot(timeInstants, median(speeds{i},2,'omitnan'),'b', color=[0.5,0.5,1]);
end
l1=plot(timeInstants, median(speeds{1},2,'omitnan'),'b',LineWidth=2);
l2=plot(timeInstants, median(speed_sim,2,'omitnan'),'k',LineWidth=2);
xlim([0,max(timeInstants)])
xlabel('$t$ [s]','Interpreter','Latex','FontSize',16)
ylabel('$v$ [px/s]','Interpreter','Latex','FontSize',16)
rng=ylim;
ylim([0,100])
if isvarname('u')
    highlightInputs(timeInstants, u, 'r', 0.25)
end
legend([l1,l2],'REAL','SIMULATED')
box on
subplot(2,1,2)
hold on
for i=2:length(data_folders)
    plot(timeInstants(1:end-1), median(abs(omegas{i}),2,'omitnan'),'b', color=[0.5,0.5,1]);
end
l1=plot(timeInstants(1:end-1), median(abs(omegas{1}),2,'omitnan'),'b',LineWidth=2);
l2=plot(timeInstants(1:end-1), median(abs(omega_sim),2,'omitnan'),'k',LineWidth=2);
xlabel('$t$ [s]','Interpreter','Latex','FontSize',16)
ylabel('$|\omega|$ [rad/s]','Interpreter','Latex','FontSize',16)
rng=ylim;
ylim([0,1.5])
if isvarname('u')
    highlightInputs(timeInstants, u, 'r', 0.25)
end
legend([l1,l2],'REAL','SIMULATED')
box on
if outputDir
    saveas(gcf,fullfile(data_folders(1), 'comparison_time_plot'))
    saveas(gcf,fullfile(data_folders(1), 'comparison_time_plot'),'png')
end

figure % NMSE
hold on
rng = max([nmse_speed_med,nmse_omega_med],[],'all');
plot(ones(size(nmse_speed_med(2:end)))*1, nmse_speed_med(2:end), 'o','color', [0.5,0.5,1])
plot([1], nmse_speed_med(1), 'o','color', 'b', 'LineWidth', 2)
plot(ones(size(nmse_omega_med(2:end)))*2, nmse_omega_med(2:end), 'o','color', [1,0.5,0.5])
plot([2], nmse_omega_med(1), 'o','color', 'r', 'LineWidth', 2)
plot(ones(size(nmse_total_med(2:end)))*3, nmse_total_med(2:end), 'o','color', [0.5,0.5,0.5])
plot([3], nmse_total_med(1), 'o','color', 'k', 'LineWidth', 2)
title('NMSE')
xlim([0,4])



