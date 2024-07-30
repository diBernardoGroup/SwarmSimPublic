clear
close all

% Add the path to your data
% experiments_folder="C:\Users\david\OneDrive - Università di Napoli Federico II\Research\Data\DOME\";    % DAVIDE
experiments_folder="/Volumes/DOMEPEN/Experiments/comparisons";                                          % ANDREA

experiment_name=[fullfile("Euglena_switch_10","combo5")];
% experiment_name=[fullfile("Volvox_switch_10","combo5")];

identification_file_name = 'identification_GB_absw_meaninit.txt';
identification_method = 'OLS+GB'; %OLS+GB
downSampling = 1;

% EUGLENA PARAMETERS
min_duration = 10; %[s]
no_mu_w = true;
use_wabs = true;   % use abs(w) for identification, then recover parameters for w
%parameters [theta, alpha, beta, mu]
init_v    = 'identification_GB_mean.txt';
init_wabs = 'identification_GB_mean.txt';
% init_v    = []; %[0.15, 0, -45, 75];
% init_wabs = []; %[0.15, 0, 0.6,  0];
% limits_v = [init_v;init_v]; 
% limits_w = [init_w;init_w];
limits_v =    [  0    0 -inf   0; 
               inf    0   0   inf]; 
limits_wabs = [  0    0   0    0; 
               inf    0  inf  inf];
           
% % VOLVOX PARAMETERS
% min_duration = 10; %[s]
% no_mu_w = false;
% use_wabs = false;  % use abs(w) for identification, then recover parameters for w
% %parameters [theta, alpha, beta, mu]
% % init_v    = 'identification_GB_mean.txt';
% % init_w    = 'identification_GB_mean.txt';
% init_v    = []; %[0.15, 0, -45, 50];
% init_w    = []; %[0.15, 0,   0,  0];
% limits_v =    [  0    0 -inf   0; 
%                inf    0   0   inf]; 
% limits_w =    [  0    0   0    0; 
%                inf    0   0   inf];
        
deltaT = 0.5;        % sampling time step
dT = 0.01;           % integration time step (for simulation)
thresholdfactor = 1; % parameters outliers detection

current_folder = fileparts(which('AnalyseDOMEexp'));
addpath(genpath(current_folder));

data_folder=fullfile(experiments_folder,experiment_name);

%% Load data
% load init values from previous identification
if exist('init_v','var') && ischar(init_v)
    fprintf('Load init values of v from %s.\n', init_v)
    identification=readtable(fullfile(data_folder,init_v));
    init_v = median([identification.theta_s, identification.alpha_s, identification.beta_s, identification.mu_s],1);
end
if exist('init_w','var') && ischar(init_w)
    fprintf('Load init values of omega from %s.\n', init_w)
    identification=readtable(fullfile(data_folder,init_w));
    init_w = median([identification.theta_w, identification.alpha_w, identification.beta_w, identification.mu_w],1);
end
if exist('init_wabs','var') && ischar(init_wabs)
    fprintf('Load init values of abs(omega) from %s.\n', init_wabs)
    identification=readtable(fullfile(data_folder,init_wabs));
    init_wabs = median([identification.theta_w, identification.alpha_w, identification.beta_w, identification.sigma_w ./ sqrt(identification.theta_w*pi)],1);
end

% LOAD DATA AND INPUTS
%identification=readtable(fullfile(current_folder,'identification.txt'));
speed  = load(fullfile(data_folder,'speeds_smooth.txt'));
omega  = load(fullfile(data_folder,'ang_vel_smooth.txt'));
inputs = load(fullfile(data_folder,'inputs.txt'));
u=inputs(:,1)/255;

% SMOOTH AND AVERAGE DATA
% speed = movmean(speed,5,'omitnan');
% omega = movmean(omega,5,'omitnan');

% speed = median(speed,2,'omitnan');
% omega = median(abs(omega),2,'omitnan');
% omega = median(omega,2,'omitnan');
% speed = mean(speed,2,'omitnan');
% omega = mean(abs(omega),2,'omitnan');
% omega = mean(omega,2,'omitnan');

% speed = movmean(speed,5,'omitnan');
% omega = movmean(omega,5,'omitnan');

% DOWNSAMPLING
speed  = speed(1:downSampling:end,:);
omega  = omega(1:downSampling:end,:);
u = u(1:downSampling:end);
deltaT = deltaT*downSampling;

N = size(speed,2);
timeInstants = [0:size(speed,1)-1] * deltaT;
agents = [0:N-1]';

% COMPUTE INPUT DERIVATIVE
u_dot_BE = [0;diff(u)]/deltaT;
u_dot_grad = gradient(u)/deltaT;
%Setting the right derivative 
u_dot = u_dot_BE;
u_dot = max(u_dot,0);
u_matrix = [u, u_dot];

u_signw={};
for i=1:N
    u_signw{i} = u_matrix(1:size(omega,1),:) .* sign(omega(:,i));
end

%% Identification
tic
if strcmp(identification_method,'GA')
    %Identification of v
    [mu_s, theta_s, sigma_s, gains_s] = SDE_parameters_est(speed, u_matrix, deltaT, identification_method, min_duration, false, @id_fcn_v, limits_v);
    %Identification of w
    [mu_w, theta_w, sigma_w, gains_w] = SDE_parameters_est(omega, u_matrix, deltaT, identification_method, min_duration, no_mu_w, @id_fcn_w, limits_w);
else
    %Identification of v
    fprintf('Identifying v ')
    [mu_s, theta_s, sigma_s, gains_s] = SDE_parameters_est(speed, u_matrix, deltaT, identification_method, min_duration, false, 'continouosSys', limits_v, init_v);
    
    %Identification of w
    fprintf('Identifying omega ')
    if use_wabs
        [mu_wabs, theta_w, sigma_w, gains_w] = SDE_parameters_est(abs(omega), u_matrix, deltaT, identification_method, min_duration, false, 'continouosSys', limits_wabs, init_wabs); % abs omega
        sigma_w = mu_wabs .* sqrt(theta_w*pi); % from std of abs(w) to std of w
        mu_w = mu_wabs * 0;                    % from mean of abs(w) to mean of w
    else
        %[mu_w, theta_w, sigma_w, gains_w] = SDE_parameters_est(omega, u_matrix, deltaT, identification_method, min_duration, no_mu_w, 'continouosSys_sgn', limits_w);
        [mu_w, theta_w, sigma_w, gains_w] = SDE_parameters_est(omega, u_signw,  deltaT, identification_method, min_duration, no_mu_w, 'continouosSys', limits_w, init_w);
    end
end
toc

%% Save data
%Approximate to the 4th figure
mu_s=round(mu_s,4);
theta_s=round(theta_s,4);
sigma_s=round(sigma_s,4);
mu_w=round(mu_w,4);
theta_w=round(theta_w,4);
sigma_w=round(sigma_w,4);
gains_w=round(gains_w,4); 
gains_s=round(gains_s,4);

% Gains(1)=alpha Gains(2)=beta
alpha_s = gains_s(:,1); 
beta_s = gains_s(:,2);
alpha_w = gains_w(:,1);
beta_w = gains_w(:,2);

%Compute average and standard deviation of speed and angular velocity
mean_s  = round(mean(speed,'omitnan')',4);
std_s   = round( std(speed,'omitnan')',4);
mean_w  = round(mean(omega,'omitnan')',4);
std_w   = round( std(omega,'omitnan')',4);

%Save the identification data 
identification = table(agents, mu_s, theta_s, sigma_s, alpha_s, beta_s, mu_w, theta_w, sigma_w, alpha_w, beta_w, mean_s, std_s, mean_w, std_w);
nan_ids = isnan(identification.mu_s) | isnan(identification.mu_w);
identification(nan_ids,:) = [];

% must be positive
for i=["mu_s","theta_s","sigma_s","theta_w","sigma_w"] %,"mu_w"
    identification(identification.(i) <= 0,:) = [];
end
% reject outliers
for i=["mu_s","theta_s","sigma_s","alpha_s","beta_s","mu_w","theta_w","sigma_w","alpha_w","beta_w"]
    identification(isoutlier(identification.(i),'quartiles',thresholdfactor=thresholdfactor),:) = [];
end
disp(['identified ',num2str(size(identification,1)),' valid agents out of ',num2str(length(agents))])
disp([num2str(sum(isnan(mean(speed,'omitnan')))),' empty tracks, ', num2str(sum(nan_ids)-sum(isnan(mean(speed,'omitnan')))),' removed during identification, ' num2str(length(agents)-size(identification,1)-sum(nan_ids)),' removed because non valid values or outliers'])


writetable(identification,fullfile(data_folder, identification_file_name) ,'Delimiter',' ')

%% simulate average behaviour
t_sim=0:dT:max(timeInstants);
s_sim=nan(length(t_sim),1);
w_sim=nan(length(t_sim),1);
u_sim=nan(length(t_sim),2);
s_sim(1)=median(identification.mu_s);
w_sim(1)=median(abs(omega),'all','omitnan');
theta_s_med = median(identification.theta_s);
mu_s_med    = median(identification.mu_s);
alpha_s_med = median(identification.alpha_s);
beta_s_med  = median(identification.beta_s);
theta_w_med = median(identification.theta_w);
mu_w_med    = median(identification.sigma_w ./ sqrt(identification.theta_w*pi));
alpha_w_med = median(identification.alpha_w);
beta_w_med  = median(identification.beta_w);
for i=1:length(t_sim)-1
    u_sim(i,:)= [u(ceil(i*dT/deltaT)),u_dot(ceil(i*dT/deltaT))]; 
%     u_sim(i,:)= [interp1(timeInstants,u,i*dT),interp1(timeInstants,u_dot,i*dT)]; 
%     u_sim(i,1)= interp1(timeInstants,u,i*dT);
%     if i==1
%         u_sim(i,2)=0; 
%     else
%         u_sim(i,2)=(u_sim(i,1)-u_sim(i-1,1))/dT;
%     end
    s_sim(i+1)= s_sim(i) + (theta_s_med*(mu_s_med-s_sim(i)) + alpha_s_med*u_sim(i,1) + beta_s_med*u_sim(i,2) )*dT;
    w_sim(i+1)= w_sim(i) + (theta_w_med*(mu_w_med-w_sim(i)) + alpha_w_med*u_sim(i,1) + beta_w_med*u_sim(i,2) )*dT;
end

mse_speed = mean((mean(speed,2,'omitnan')-interp1(t_sim,s_sim,timeInstants)').^2);
mse_omega = mean((mean(abs(omega),2,'omitnan')-interp1(t_sim,w_sim,timeInstants)').^2);
nmse_speed = goodnessOfFit(interp1(t_sim,s_sim,timeInstants)', mean(speed,2,'omitnan'), 'NMSE');
nmse_omega = goodnessOfFit(interp1(t_sim,w_sim,timeInstants)', mean(abs(omega),2,'omitnan'), 'NMSE');
nmse_total = mean([nmse_speed, nmse_omega]);
%disp(['MSE from mean for speed: ',num2str(mse_speed),' and omega: ',num2str(mse_omega)])
disp(['NMSE from mean for speed:',num2str(nmse_speed),' and omega: ',num2str(nmse_omega),' total: ', num2str(nmse_total)])

mse_speed = mean((median(speed,2,'omitnan')-interp1(t_sim,s_sim,timeInstants)').^2);
mse_omega = mean((median(abs(omega),2,'omitnan')-interp1(t_sim,w_sim,timeInstants)').^2);
nmse_speed = goodnessOfFit(interp1(t_sim,s_sim,timeInstants)', median(speed,2,'omitnan'), 'NMSE');
nmse_omega = goodnessOfFit(interp1(t_sim,w_sim,timeInstants)', median(abs(omega),2,'omitnan'), 'NMSE');
nmse_total = mean([nmse_speed, nmse_omega]);
disp(['NMSE from median for speed:',num2str(nmse_speed),' and omega: ',num2str(nmse_omega),' total: ', num2str(nmse_total)])
disp(['Saved as ',identification_file_name])

%% PLOTS
[~,f,~]=fileparts(identification_file_name);
if ~exist(fullfile(data_folder,'plots'),'dir'); mkdir(fullfile(data_folder,'plots')); end

figure % PARAMETERS BOXPLOTS
for i=1:10 
    ax=subplot(2,5,i);
    myboxplot(identification(:,i+1), false, 3);
    set(ax,'PositionConstraint','innerposition')
    yline(0,'Color',[0.5,0.5,0.5])
    %l=max(identification.(i+1))*1.1;ylim([-l/15,l]);yticks([0:l/3:l])
    set(gca,'FontSize',16)
end
saveas(gcf,fullfile(data_folder,'plots',[f '_parameters']))
saveas(gcf,fullfile(data_folder,'plots',[f '_parameters']),'png')

figure % PARAMETERS HISTOGRAMS
for i=1:10 
    ax=subplot(2,5,i);
    histogram(identification{:,i+1});
    xlabel(identification.Properties.VariableNames{i+1})
    set(ax,'PositionConstraint','innerposition')
    yline(0,'Color',[0.5,0.5,0.5])
    %l=max(identification.(i+1))*1.1;ylim([-l/15,l]);yticks([0:l/3:l])
    set(gca,'FontSize',16)
end
saveas(gcf,fullfile(data_folder,'plots',[f '_parameters_hist']))
saveas(gcf,fullfile(data_folder,'plots',[f '_parameters_hist']),'png')

figure % TIME PLOT - SPEED and ANGULAR VELOCITY
subplot(2,4,[1 2 3])
line1=plotWithShade(timeInstants, median(speed,2,'omitnan'), min(speed, [], 2,'omitnan'), max(speed, [], 2,'omitnan'), 'b', 0.3);
line2=plot(timeInstants, mean(speed,2,'omitnan'),'g',LineWidth = 2);
line3=plot(t_sim,s_sim,'r',LineWidth = 1);
highlightInputs(timeInstants, u, 'r', 0.25)
legend([line1,line2,line3],'experiment (median)', 'experiment (mean)', 'simulation')
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
line1=plotWithShade(timeInstants, median(abs(omega),2,'omitnan'), min(abs(omega), [], 2,'omitnan'), max(abs(omega), [], 2,'omitnan'), 'b', 0.3);
line2=plot(timeInstants, mean(abs(omega),2,'omitnan'),'g',LineWidth = 2);
line3=plot(t_sim,abs(w_sim),'r',LineWidth = 1);
highlightInputs(timeInstants, u, 'r', 0.25)
legend([line1,line2,line3],'experiment (median)', 'experiment (mean)', 'simulation')
xlabel('t [s]')
ylabel('ang. vel. [rad/s]')
xlim([0,inf])
rng=ylim;
box on
subplot(2,4,8)
h=histogram(abs(omega(:)),'Orientation','horizontal');
ylim(rng);
set(gca,'xtick',[])
saveas(gcf,fullfile(data_folder,'plots',[f '_sim']))
saveas(gcf,fullfile(data_folder,'plots',[f '_sim']),'png')

% figure
% corrplot(identification(:,2:11))
% saveas(gcf,fullfile(data_folder,'plots',[f '_corrplot']))
% saveas(gcf,fullfile(data_folder,'plots',[f '_corrplot']),'png')
% 
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
