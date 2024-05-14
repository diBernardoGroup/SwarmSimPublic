clear
close all

% experiments data folders
experiments_folder = '/Volumes/DOMEPEN/Experiments';

% experiments must have the same inputs
% switch 10s all
tag = 'Euglena_switch_10';
experiments_names = ["2023_06_15_Euglena_7", "2023_06_26_Euglena_23","2023_06_26_Euglena_24", "2023_07_10_Euglena_15", "2023_07_10_Euglena_16"];

% switch 10s selected
tag = 'Euglena_switch_10';
experiments_names = ["2023_06_15_Euglena_7", "2023_06_26_Euglena_24", "2023_07_10_Euglena_15"];

% switch 5s
tag = 'Euglena_switch_5';
experiments_names = ["2023_06_15_Euglena_8", "2023_06_26_Euglena_25", "2023_07_10_Euglena_18"];

% switch 1s
tag = 'Euglena_switch_1';
experiments_names = ["2023_06_15_Euglena_11", "2023_06_26_Euglena_28", "2023_07_10_Euglena_19"];

% off
tag = 'Euglena_OFF';
experiments_names = ["2023_06_15_Euglena_1", "2023_06_26_Euglena_13", "2023_07_10_Euglena_6"];

% OFF-ON-OFF 75
tag = 'Euglena_75_ON';
experiments_names = ["2023_06_15_Euglena_2", "2023_06_26_Euglena_16", "2023_07_10_Euglena_8"]; % old
experiments_names = ["2023_06_15_Euglena_2", "2023_06_26_Euglena_15", "2023_07_10_Euglena_8"]; % old

% % OFF-ON-OFF 150
% tag = 'Euglena_150_ON';
% experiments_names = ["2023_06_15_Euglena_3", "2023_06_26_Euglena_18", "2023_07_10_Euglena_10"];
% 
% % OFF-ON-OFF 255
% tag = 'Euglena_255_ON';
% experiments_names = ["2023_06_15_Euglena_4", "2023_06_26_Euglena_20", "2023_07_10_Euglena_12"];
% 
% % Ramp
% tag = 'Euglena_ramp';
% experiments_names = ["2023_06_15_Euglena_6", "2023_06_26_Euglena_22", "2023_07_10_Euglena_14"];

output_folder = fullfile("/Volumes/DOMEPEN/Experiments/comparisons/",tag,"/combo");

%% LOAD AND CONCATENATE DATA
speed = [];
omega = [];

for i=1:length(experiments_names)
    tracking = getLastTracking(fullfile(experiments_folder,experiments_names(i)));
    data_folders(i) =  fullfile(experiments_folder,experiments_names(i), tracking);
    speed  = horzcat(speed, load(fullfile(data_folders(i),'speeds_smooth.txt')));
    omega  = horzcat(omega, load(fullfile(data_folders(i),'ang_vel_smooth.txt')));
    
    if i==1
        inputs = load(fullfile(data_folders(i),'inputs.txt'));
    else
        assert( all(inputs==load(fullfile(data_folders(i),'inputs.txt')),'all'), 'Experiments have different inputs!')
    end
end

%% SAVE DATA
if ~exist(output_folder,'dir'); mkdir(output_folder); end
if ~exist(fullfile(output_folder,'plots'),'dir'); mkdir(fullfile(output_folder,'plots')); end

disp(['Data saved in ',char(output_folder)])

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
xlim([0,max(timeInstants)])
ylim([0,120])
if isvarname('u')
    highlightInputs(timeInstants, u, 'r', 0.25)
end
l1=plotWithShade(timeInstants, median(speed,2,'omitnan'), quantile(speed,0.1,2), quantile(speed,0.9,2), 'b', 0.3);
xlabel('$t$ [s]','Interpreter','Latex','FontSize',16)
ylabel('$v$ [$\mu$m/s]','Interpreter','Latex','FontSize',16)
%legend([l1,l2],'REAL','SIMULATED')
rng=ylim;
box on
subplot(2,1,2)
xlim([0,max(timeInstants)])
ylim([0,2])
if isvarname('u')
    highlightInputs(timeInstants, u, 'r', 0.25)
end
l1=plotWithShade(timeInstants, median(abs(omega),2,'omitnan'), quantile(abs(omega),0.1,2), quantile(abs(omega),0.9,2), 'b', 0.3);
xlabel('$t$ [s]','Interpreter','Latex','FontSize',16)
ylabel('$|\omega|$ [rad/s]','Interpreter','Latex','FontSize',16)
%legend([l1,l2],'REAL','SIMULATED')
rng=ylim;
box on
saveas(gcf,fullfile(output_folder, 'plots', 'time_plot'))
saveas(gcf,fullfile(output_folder, 'plots', 'time_plot'),'png')


