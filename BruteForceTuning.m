%
%BruteForceTuning Set the parameters and ranges of the control gains.
%   Multiple simulations are executed to find the best values of the control gains.
%
%   Notes: 
%       Running this script can take long time (up to days)
%       For better performances install the parallel computing toolbox
%       This script can be modified to vary any two parameters, besides the
%       control gains, such as the parameters of the interaction functions
%
%   See also: Launcher, SequentialLauncher
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%

%% Clear environment
clear
close all
clc

%% Parameters

Ntimes=4;       % How many times simulator is launched on each config (couple of gains)

N=100;          %number of agents (N)
LinkNumber=4;   %number of links (6=triangular lattice, 4=square lattice, 3=hexagonal lattice) (L)

% descrption of the radial interaction function
IntFunctionStruct=struct('function','Lennard-Jones','parameters',[0.15, 5]);

% adaptation gains
alpha = 0;
beta = 0;

% thresholds
regularity_thresh=0.2;      % threshold value for regularity metrics (e^*_theta)
compactness_thresh=0.3;     % threshold value for compactness metrics (e^*_L)
length_thresh=0.1;          % threshold value for error in the lenght of the links

Tmax=200;    % maximum simulation time (simulation is stopped earlier if steady state is reached)

sigma = 0;   % standard deviation of noise

MaxSensingRadius=inf;   % sensing radius of the agents (R_s)

%output options
drawON=false;       % draw swarm during simulation (set to false for extensive simulations)
getMetrics=true;    % acquire metrics during the simulation (getMetrics=false discard settling times and stop times)

% robustness tests
AgentsRemoval=false;
NoiseTest=false;
dynamicLattice=false;

%ranges of values of control gains
G_r_min = 0;
G_r_max = 30;
G_n_min = 0;
G_n_max = 30;

%grid of values
step=1;
G_r_vec = [G_r_min:step:G_r_max];
G_n_vec = [G_n_min:step:G_n_max];

% values to import from each simulation
e_theta = zeros(1, Ntimes);
e_L = zeros(1, Ntimes);
e_d = zeros(1, Ntimes);
Tr_vec = zeros(1, Ntimes);
success_vec = zeros(1, Ntimes);
stopTime_vec = zeros(1, Ntimes);

e_theta_mean = zeros(length(G_r_vec), length(G_n_vec));
e_L_mean = zeros(length(G_r_vec), length(G_n_vec));
e_d_mean = zeros(length(G_r_vec), length(G_n_vec));
Tr_mean = zeros(length(G_r_vec), length(G_n_vec));
success_mean = zeros(length(G_r_vec), length(G_n_vec));
stopTime_mean = zeros(length(G_r_vec), length(G_n_vec));
stopTime_max = zeros(length(G_r_vec), length(G_n_vec));

%% Create initial conditions
rng(0,'twister'); % reproducible results
x0Data=nan(Ntimes,N,2);
for k_times=1:Ntimes
    x0Data(k_times,:,:)=randCircle(N, 2); % circle of radius 2
    %x0Data(k_times,:,:)=randCircle(N, mod(k_times*2-1,8)+1); % circles of radius 2,4,6 ed 8
end

%% Simulate Ntimes over all values of control gains
for i_times= 1:length(G_r_vec)
    G_r = G_r_vec(i_times); 
    for j_times= 1:length(G_n_vec)
        G_n = G_n_vec(j_times);
        disp('Simulation ' + string((i_times-1)*length(G_n_vec)+j_times) + ' of ' + string(length(G_r_vec)*length(G_n_vec)))

        tic
        parfor k_times=  1:Ntimes
            x0=squeeze(x0Data(k_times,:,:));
                
            %Gr-Gn tuning
            [T_r, success, final_e_theta, final_e_L, final_e_d, finalGRadial, finalGNormal, stopTime] = Simulator(x0, LinkNumber, G_r, G_n, regularity_thresh, compactness_thresh, Tmax, sigma, drawON, getMetrics, IntFunctionStruct, AgentsRemoval, NoiseTest, MaxSensingRadius, alpha, beta, dynamicLattice);
            
            %import data from the simulation
            e_theta(k_times) = final_e_theta;
            e_L(k_times) = final_e_L;
            e_d(k_times) = final_e_d;
            Tr_vec (k_times) = T_r;
            success_vec(k_times) = success;
            stopTime_vec(k_times) = stopTime;
        end
        toc
        
        e_theta_mean(i_times, j_times) = mean(e_theta);
        e_L_mean(i_times, j_times) = mean(e_L);
        e_d_mean(i_times, j_times) = mean(e_d);
        Tr_mean(i_times, j_times) = mean(Tr_vec, 'omitnan');
        success_mean(i_times, j_times) = mean(success_vec);
        stopTime_mean(i_times, j_times) = mean(stopTime_vec, 'omitnan');
        stopTime_max(i_times, j_times) = max(stopTime_vec,[], 'includenan');
    end
end

window=1;
filtered_Tr=movmean2(Tr_mean,window);
filtered_e_theta=movmean2(e_theta_mean,window);
filtered_e_L=movmean2(e_L_mean,window);
filtered_e_d=movmean2(e_d_mean,window);
filtered_succ=movmean2(success_mean,window);

% find optimal gains
[row, col, costMap]=findOptim(cat(3, filtered_e_theta/regularity_thresh, filtered_e_L/compactness_thresh));
optimalGains=[G_r_vec(row),G_n_vec(col)];

%% Plots

%e_theta
figure
[~,lplot]=mysurfc(G_r_vec, G_n_vec, filtered_e_theta, regularity_thresh, 0.6);
xlabel('G_r')
ylabel('G_n')
title('$\bar{e}_{\theta}$', 'interpreter', 'latex')
hold on
optimplot=scatter3(G_r_vec(row),G_n_vec(col),10*ones(length(col),1), 100,'black','filled');
xlim([-inf, inf])
ylim([-inf, inf])
set(gca, 'XTick', sort(unique([G_r_min, G_r_max, get(gca, 'XTick')])));
set(gca, 'YTick', sort(unique([G_n_min, G_n_max, get(gca, 'YTick')])));
set(gca,'FontSize',14)
legend([lplot,optimplot],{'$e_{\theta}^*=$'+string(regularity_thresh),'Optimal Gains ('+string(optimalGains(1))+', '+string(optimalGains(2))+')'}, 'Interpreter', 'latex','FontSize',14)

%e_L
figure
[~,lplot]=mysurfc(G_r_vec, G_n_vec,filtered_e_L,compactness_thresh, 0.15);
xlabel('G_r')
ylabel('G_n')
title('$\bar{e}_{L}$', 'interpreter', 'latex')
hold on
optimplot=scatter3(G_r_vec(row),G_n_vec(col),10*ones(length(col),1), 100,'black','filled');
xlim([-inf, inf])
ylim([-inf, inf])
set(gca, 'XTick', sort(unique([G_r_min, G_r_max, get(gca, 'XTick')])));
set(gca, 'YTick', sort(unique([G_n_min, G_n_max, get(gca, 'YTick')])));
set(gca,'FontSize',14)
legend([lplot,optimplot],{'$e_{L}^*=$'+string(compactness_thresh),'Optimal Gains ('+string(optimalGains(1))+', '+string(optimalGains(2))+')'}, 'Interpreter', 'latex','FontSize',14)

% %Tr
% figure
% [~,lplot]=mysurfc(G_r_vec, G_n_vec, filtered_Tr, 20, 0.15);
% xlabel('G_r')
% ylabel('G_n')
% title('$T_{r}$', 'interpreter', 'latex')
% hold on
% optimplot=scatter3(G_r_vec(row),G_n_vec(col),10*ones(length(col),1), 100,'black','filled');
% xlim([-inf, inf])
% ylim([-inf, inf])
% set(gca, 'XTick', sort(unique([G_r_min, G_r_max, get(gca, 'XTick')])));
% set(gca, 'YTick', sort(unique([G_n_min, G_n_max, get(gca, 'YTick')])));
% set(gca,'FontSize',14)
% legend([lplot,optimplot],{'$20$','Optimal Gains ('+string(optimalGains(1))+', '+string(optimalGains(2))+')'}, 'Interpreter', 'latex','FontSize',14)

% Cost Map
figure
[~,lplot]=mysurfc(G_r_vec, G_n_vec,costMap,1, 0.15);
caxis([0 2])                 
xlabel('G_r')
ylabel('G_n')
title('Cost')
hold on
optimplot=scatter3(optimalGains(1),optimalGains(2),10*ones(length(col),1), 100,'black','filled');
xlim([-inf, inf])
ylim([-inf, inf])
set(gca, 'XTick', sort(unique([G_r_min, G_r_max, get(gca, 'XTick')])));
set(gca, 'YTick', sort(unique([G_n_min, G_n_max, get(gca, 'YTick')])));
set(gca,'FontSize',14)
legend([lplot,optimplot],{'$1$','Optimal Gains ('+string(optimalGains(1))+', '+string(optimalGains(2))+')'}, 'Interpreter', 'latex','FontSize',14)

% Success Rate Map
figure
[~,lplot]=mysurfc(G_r_vec, G_n_vec,filtered_succ,0.75, 0.15);
caxis([0 1])                 
xlabel('G_r')
ylabel('G_n')
title('Success Rate')
hold on
optimplot=scatter3(G_r_vec(row),G_n_vec(col),10*ones(length(col),1), 100,'black','filled');
xlim([-inf, inf])
ylim([-inf, inf])
set(gca, 'XTick', sort(unique([G_r_min, G_r_max, get(gca, 'XTick')])));
set(gca, 'YTick', sort(unique([G_n_min, G_n_max, get(gca, 'YTick')])));
set(gca,'FontSize',14)
legend([lplot,optimplot],{'$0.75$','Optimal Gains ('+string(optimalGains(1))+', '+string(optimalGains(2))+')'}, 'Interpreter', 'latex','FontSize',14)

% % Stop Time Map
% figure
% heatmap( G_r_vec, fliplr(G_n_vec), fliplr(stopTime_max)')
% caxis([0 Tmax*0.5])                 
% xlabel('G_r')
% ylabel('G_n')
% title('Stop Time Max')
% set(gca,'FontSize',14)
% 
% % Stop Time Map
% figure
% heatmap( G_r_vec, fliplr(G_n_vec), fliplr(stopTime_mean)')
% caxis([0 Tmax*0.25])                 
% xlabel('G_r')
% ylabel('G_n')
% title('Stop Time Mean')
% set(gca,'FontSize',14)

% %e_d
% figure
% [~,lplot]=mysurfc(G_r_vec, G_n_vec,filtered_e_d,length_thresh, 0.15);
% xlabel('G_r')
% ylabel('G_n')
% title('$\bar{e}_{d}$', 'interpreter', 'latex')
% hold on
% optimplot=scatter3(G_r_vec(row),G_n_vec(col),10*ones(length(col),1), 100,'black','filled');
% xlim([-inf, inf])
% ylim([-inf, inf])
% set(gca, 'XTick', sort(unique([G_r_min, G_r_max, get(gca, 'XTick')])));
% set(gca, 'YTick', sort(unique([G_n_min, G_n_max, get(gca, 'YTick')])));
% set(gca,'FontSize',14)
% legend([lplot,optimplot],{'$e_{d}^*=$'+string(length_thresh),'Optimal Gains ('+string(optimalGains(1))+', '+string(optimalGains(2))+')'}, 'Interpreter', 'latex','FontSize',14)

