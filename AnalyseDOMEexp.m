clear
close all


current_folder = '/Volumes/DOMEPEN/Experiments/2023_06_15_Euglena_7/tracking_2023_10_16';

deltaT = 0.5;
thresholdfactor = 5;

%% Load data
%identification=readtable(fullfile(current_folder,'identification.txt'));
speed=load(fullfile(current_folder,'speeds_smooth.txt'));
omega=load(fullfile(current_folder,'ang_vel_smooth.txt'));
inputs=load(fullfile(current_folder,'inputs.txt'));

N = size(speed,2);
agents = [0:N-1]';
u=inputs(:,1)/255;

%% Identification
[mu_s, theta_s, sigma_s, alpha_s] = SDE_parameters_est(speed, u, deltaT, 'OLS');
[mu_w, theta_w, sigma_w, alpha_w] = SDE_parameters_est(omega, u, deltaT,'OLS');
mu_s=round(mu_s,4); theta_s=round(theta_s,4); sigma_s=round(sigma_s,4); alpha_s=round(alpha_s,4); mu_w=round(mu_w,4); theta_w=round(theta_w,4); sigma_w=round(sigma_w,4); alpha_w=round(alpha_w,4);  

mean_s  = round(mean(speed,'omitnan')',4);
std_s   = round(std(speed,'omitnan')',4);
mean_w  = round(mean(omega,'omitnan')',4);
std_w   = round(std(omega,'omitnan')',4);

identification = table(agents, mu_s, theta_s, sigma_s, alpha_s, mu_w, theta_w, sigma_w, alpha_w, mean_s, std_s, mean_w, std_w);
nan_ids = isnan(identification.mu_s) | isnan(identification.mu_w);
identification(nan_ids,:) = [];
for i=["mu_s","theta_s","sigma_s","theta_w","sigma_w"]
    identification(identification.(i) < 0,:) = [];
end
for i=["mu_s","theta_s","sigma_s","alpha_s","theta_w","sigma_w","alpha_w"]
    identification(isoutlier(identification.(i),'quartiles',thresholdfactor=thresholdfactor),:) = [];
end
disp(['identified ',num2str(size(identification,1)),' valid agents out of ',num2str(length(agents))])
disp([num2str(sum(nan_ids)),' removed because nans, ' num2str(length(agents)-size(identification,1)-sum(nan_ids)),' removed because negative or outliers'])


writetable(identification,fullfile(current_folder, 'identification.txt') ,'Delimiter',' ')

%% PLOTS
timeInstants = [1:size(speed,1)] * deltaT;

figure
for i=1:8 
subplot(1,8,i)
myboxplot(identification(:,i+1), false, 3);
yline(0,'Color',[0.5,0.5,0.5])
%l=max(identification.(i+1))*1.1;ylim([-l/15,l]);yticks([0:l/3:l])
set(gca,'FontSize',16)
end

figure % TIME PLOT - SPEED and ANGULAR VELOCITY
subplot(2,4,[1 2 3])
plotWithShade(timeInstants, median(speed,2,'omitnan'), min(speed, [], 2,'omitnan'), max(speed, [], 2,'omitnan'), 'b', 0.3);
highlightInputs(timeInstants, u, 'r', 0.25)
xlabel('t [s]')
ylabel('speed')
xlim([0,inf])
rng=ylim;
box on
subplot(2,4,4)
h=histogram(speed(:),'Orientation','horizontal');
ylim(rng);
set(gca,'xtick',[])
subplot(2,4,[5 6 7])
plotWithShade(timeInstants, median(abs(omega),2,'omitnan'), min(abs(omega), [], 2,'omitnan'), max(abs(omega), [], 2,'omitnan'), 'b', 0.3);
highlightInputs(timeInstants, u, 'r', 0.25)
xlabel('t [s]')
ylabel('ang. vel. [rad/s]')
rng=ylim;
box on
subplot(2,4,8)
h=histogram(abs(omega(:)),'Orientation','horizontal');
ylim(rng);
set(gca,'xtick',[])