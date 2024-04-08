clear
close all

% experiments must have the same inputs
data_folders = ["/Volumes/DOMEPEN/Experiments/2023_06_15_Euglena_7/tracking_2023_10_16",
    "/Volumes/DOMEPEN/Experiments/2023_06_26_Euglena_23/tracking_2023_10_12",
    "/Volumes/DOMEPEN/Experiments/2023_06_26_Euglena_24/tracking_2024_04_08",
    "/Volumes/DOMEPEN/Experiments/2023_07_10_Euglena_15/tracking_2023_10_12",
    "/Volumes/DOMEPEN/Experiments/2023_07_10_Euglena_16/tracking_2024_04_08"]; % experiments data folders

output_folder = "/Volumes/DOMEPEN/Experiments/comparisons/Euglena_switch_10/combo";

%% LOAD AND CONCATENATE DATA
speed = [];
omega = [];

for i=1:length(data_folders)
    speed  = horzcat(speed, load(fullfile(data_folders(i),'speeds_smooth.txt')));
    omega  = horzcat(omega, load(fullfile(data_folders(i),'ang_vel_smooth.txt')));
    
    if i==1
        inputs = load(fullfile(data_folders(i),'inputs.txt'));
    else
        assert( all(inputs==load(fullfile(data_folders(i),'inputs.txt')),'all') )
    end
end

%% SAVE DATA
if ~exist(output_folder,'dir'); mkdir(output_folder); end
if ~exist(fullfile(output_folder,'plots'),'dir'); mkdir(fullfile(output_folder,'plots')); end

fileID = fopen(fullfile(output_folder, 'combo_info.txt'),'wt');
fprintf(fileID,'Combined DOME exp data\n\n');
fprintf(fileID,'Date: %s\n',datestr(now, 'dd/mm/yy'));
fprintf(fileID,'Time: %s\n\n',datestr(now, 'HH:MM'));
fprintf(fileID,'Experiments:\n');
fprintf(fileID,'%s\n',data_folders);
fclose(fileID);

writematrix(speed, fullfile(output_folder,'speeds_smooth.txt'),'Delimiter','tab');
writematrix(omega, fullfile(output_folder,'ang_vel_smooth.txt'),'Delimiter','tab');
writematrix(inputs, fullfile(output_folder,'inputs.txt'),'Delimiter','tab');

assert(all(size(load(fullfile(output_folder,'speeds_smooth.txt')))==size(speed)))

%% PLOTS
timeInstants = [0:0.5:180];
u=inputs(:,1)/255;              %select blue channel and scale in [0,1]

figure % TIME PLOT - SPEED and ANGULAR VELOCITY
subplot(2,1,1)
l1=plotWithShade(timeInstants, median(speed,2,'omitnan'), min(speed, [], 2,'omitnan'), max(speed, [], 2,'omitnan'), 'b', 0.3);
xlim([0,max(timeInstants)])
if isvarname('u')
    highlightInputs(timeInstants, u, 'r', 0.25)
end
xlabel('$t$ [s]','Interpreter','Latex','FontSize',16)
ylabel('$v$ [px/s]','Interpreter','Latex','FontSize',16)
%legend([l1,l2],'REAL','SIMULATED')
rng=ylim;
ylim([0,120])
box on
% subplot(2,4,4)
% hold on
% h=histogram(experiments{1}.speed(:),'Orientation','horizontal');
% h=histogram(experiments{2}.speed(:),'Orientation','horizontal');
% ylim(rng);
% set(gca,'xtick',[])
subplot(2,1,2)
l1=plotWithShade(timeInstants, median(abs(omega),2,'omitnan'), min(abs(omega), [], 2,'omitnan'), max(abs(omega), [], 2,'omitnan'), 'b', 0.3);
xlim([0,max(timeInstants)])
if isvarname('u')
    highlightInputs(timeInstants, u, 'r', 0.25)
end
xlabel('$t$ [s]','Interpreter','Latex','FontSize',16)
ylabel('$|\omega|$ [rad/s]','Interpreter','Latex','FontSize',16)
%legend([l1,l2],'REAL','SIMULATED')
rng=ylim;
ylim([0,3])
box on
% subplot(2,4,8)
% hold on
% h=histogram(abs(experiments{1}.omega(:)),'Orientation','horizontal', 'FaceColor', 'b', 'EdgeColor', 'none');
% h=histogram(abs(experiments{2}.omega(:)),'Orientation','horizontal', 'FaceColor', 'k', 'EdgeColor', 'none');
% ylim(rng);
% set(gca,'xtick',[])
saveas(gcf,fullfile(output_folder, 'plots', 'time_plot'))
saveas(gcf,fullfile(output_folder, 'plots', 'time_plot'),'png')


