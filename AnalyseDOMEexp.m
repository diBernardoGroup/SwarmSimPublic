clear
close all


%current_folder = '/Volumes/DOMEPEN/Experiments/2023_06_15_Euglena_1/tracking_2023_10_12'; % off
current_folder = '/Volumes/DOMEPEN/Experiments/2023_06_15_Euglena_7/tracking_2023_10_16'; % switch10s
%current_folder = '/Volumes/DOMEPEN/Experiments/2023_06_26_Euglena_19/tracking_2023_10_16'; % on255

deltaT = 0.5;
dT = 0.01;
thresholdfactor = 5;

%% Load data
%identification=readtable(fullfile(current_folder,'identification.txt'));
speed=load(fullfile(current_folder,'speeds_smooth.txt'));
omega=load(fullfile(current_folder,'ang_vel_smooth.txt'));
inputs=load(fullfile(current_folder,'inputs.txt'));

N = size(speed,2);
timeInstants = [1:size(speed,1)] * deltaT;
agents = [0:N-1]';
u=inputs(:,1)/255;
u_dot = [diff(u);0]/deltaT;

%% Identification
[mu_s, theta_s, sigma_s, gains_s] = SDE_parameters_est(speed, [u, u_dot], deltaT, 'OLS');
%[mu_w, theta_w, sigma_w, gains_w] = SDE_parameters_est(omega, [u, u_dot], deltaT,'OLS');
[mu_w, theta_w, sigma_w, gains_w] = SDE_parameters_est(abs(omega), [u, u_dot], deltaT,'OLS');
mu_s=round(mu_s,4); theta_s=round(theta_s,4); sigma_s=round(sigma_s,4); gains_s=round(gains_s,4); mu_w=round(mu_w,4); theta_w=round(theta_w,4); sigma_w=round(sigma_w,4); gains_w=round(gains_w,4);  
alpha_s = gains_s(:,1); beta_s = gains_s(:,2); alpha_w = gains_w(:,1); beta_w = gains_w(:,2);

mean_s  = round(mean(speed,'omitnan')',4);
std_s   = round( std(speed,'omitnan')',4);
mean_w  = round(mean(omega,'omitnan')',4);
std_w   = round( std(omega,'omitnan')',4);

identification = table(agents, mu_s, theta_s, sigma_s, alpha_s, beta_s, mu_w, theta_w, sigma_w, alpha_w, beta_w, mean_s, std_s, mean_w, std_w);
nan_ids = isnan(identification.mu_s) | isnan(identification.mu_w);
identification(nan_ids,:) = [];
for i=["alpha_s","beta_s","alpha_w","beta_w"]
    identification(identification.(i) == 0,:) = [];
end
for i=["mu_s","theta_s","sigma_s","theta_w","sigma_w"]
    identification(identification.(i) <= 0,:) = [];
end
for i=["mu_s","theta_s","sigma_s","alpha_s","beta_s","theta_w","sigma_w","alpha_w","beta_w"]
    identification(isoutlier(identification.(i),'quartiles',thresholdfactor=thresholdfactor),:) = [];
end
disp(['identified ',num2str(size(identification,1)),' valid agents out of ',num2str(length(agents))])
disp([num2str(sum(nan_ids)),' removed because nans, ' num2str(length(agents)-size(identification,1)-sum(nan_ids)),' removed because negative or outliers'])


writetable(identification,fullfile(current_folder, 'identification.txt') ,'Delimiter',' ')

% simulate average behaviour
s_sim=nan(max(timeInstants)/dT,1);
w_sim=nan(max(timeInstants)/dT,1);
s_sim(1)=mean(identification.mu_s);
w_sim(1)=mean(identification.mu_w);
for i=1:max(timeInstants)/dT-1
    s_sim(i+1)= s_sim(i) + (mean(identification.theta_s)*(mean(identification.mu_s)-s_sim(i)) + mean(identification.alpha_s)*u(ceil(i*dT/deltaT)) + mean(identification.beta_s)*u_dot(ceil(i*dT/deltaT)) )*dT;
    w_sim(i+1)= w_sim(i) + (mean(identification.theta_w)*(mean(identification.mu_w)-w_sim(i)) + mean(identification.alpha_w)*u(ceil(i*dT/deltaT)) + mean(identification.beta_w)*u_dot(ceil(i*dT/deltaT)) )*dT;
end

%% PLOTS

figure
for i=1:10 
    subplot(1,10,i)
    myboxplot(identification(:,i+1), false, 3);
    yline(0,'Color',[0.5,0.5,0.5])
    %l=max(identification.(i+1))*1.1;ylim([-l/15,l]);yticks([0:l/3:l])
    set(gca,'FontSize',16)
end

figure % TIME PLOT - SPEED and ANGULAR VELOCITY
subplot(2,4,[1 2 3])
plotWithShade(timeInstants, median(speed,2,'omitnan'), min(speed, [], 2,'omitnan'), max(speed, [], 2,'omitnan'), 'b', 0.3);
plot(0:dT:max(timeInstants-dT),s_sim,'r',LineWidth = 1)
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
plot(0:dT:max(timeInstants-dT),abs(w_sim),'r',LineWidth = 1)
highlightInputs(timeInstants, u, 'r', 0.25)
xlabel('t [s]')
ylabel('ang. vel. [rad/s]')
xlim([0,inf])
rng=ylim;
box on
subplot(2,4,8)
h=histogram(abs(omega(:)),'Orientation','horizontal');
ylim(rng);
set(gca,'xtick',[])