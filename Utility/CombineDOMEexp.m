clear
close all

% experiments must have the same inputs
data_folders = ["/Volumes/DOMEPEN/Experiments/2023_06_15_Euglena_7/tracking_2023_10_16",
                "/Volumes/DOMEPEN/Experiments/2023_06_26_Euglena_23/tracking_2023_10_12",
                "/Volumes/DOMEPEN/Experiments/2023_07_10_Euglena_15/tracking_2023_10_12"]; % experiments data folders

output_file = "/Volumes/DOMEPEN/Experiments/comparisons/Euglena_switch_10";

%% LOAD DATA
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

if ~exist(output_file,'dir'); mkdir(output_file); end

writematrix(speed, fullfile(output_file,'speeds_smooth.txt'),'Delimiter','tab');
writematrix(omega, fullfile(output_file,'ang_vel_smooth.txt'),'Delimiter','tab');
writematrix(inputs, fullfile(output_file,'inputs.txt'),'Delimiter','tab');

assert(all(size(load(fullfile(output_file,'speeds_smooth.txt')))==size(speed)))


