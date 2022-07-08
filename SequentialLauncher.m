%
%SequentialLauncher allows to set multiple values of the parameters and
%   launch multiple simulations for each configuration.
%   It is used to study the effect of parameters on the system.
%
%   Notes: 
%       Running this script can take long time (up to hours)
%       For better performances install the parallel computing toolbox   
%
%   See also: Launcher, BruteForceTuning
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%

%% Clear environment
clear
close all
clc

%% Parameters

Ntimes=4; % How many simulations are launched for each configuration

% default parameters

N=100;          %number of agents (N)
LinkNumber=4;   %number of links (6=triangular lattice, 4=square lattice) (L)

% descrption of the radial interaction function
IntFunctionStruct=struct('function','Lennard-Jones','parameters',[0.15, 5]);

regularity_thresh=0.2;      % threshold value for regularity metrics (e^*_theta)
compactness_thresh=0.3;     % threshold value for compactness metrics (e^*_L)

Tmax=200;    % maximum simulation time (simulation is stopped earlier if steady state is reached)

sigma = 0;   % standard deviation of noise

% control gains
G_radial= 15;   % default value for square lattice 15 (G_r)
G_normal = 8;   % default value for square lattice  8 (G_n)

% adaptation gains
alpha = 0;      
beta = 0;

MaxSensingRadius=inf;   % sensing radius of the agents (R_s)

seed=0;         % seed of the random generator (comment out to have new results on each execution)

initRadius=2;   % radius of the circle to drow the initial positions

% variable parameters
% one or multiple parameters can be modified at the same time
% all the modified parameters must have the same number of values
% the values specified here overwrite the default ones
% IntFunctionStruct cannot be modified
parameters.names={'N','initRadius'};
parameters.values={[25 50 100],[sqrt(25/25) sqrt(50/25) sqrt(100/25)]};

%% Options
% output options
drawON=false;
getMetrics=true; %getMetrics=false discard settling times and stop times

% robustness tests
AgentsRemoval = false;      % randomly remove agents during the simulation
dynamicLattice = false;     % change lattice during the simulation


%% Preallocate
Nparameters=length(parameters.names);
Nconfig=length(parameters.values{1});

e_L=zeros(1,Ntimes);    
e_theta=zeros(1,Ntimes);
Tr_vec=zeros(1,Ntimes);
success_vec=zeros(1,Ntimes);
stopTime_vec = zeros(1, Ntimes);
GNormal_vec = zeros(1, Ntimes);

e_L_mean=zeros(Nconfig,3);
e_theta_mean=zeros(Nconfig,3);
Tr_mean=zeros(1,Nconfig);
success_mean=zeros(1,Nconfig);
stopTime_max = zeros(1,Nconfig);
GNormal_mean = zeros(Nconfig,3);

%% Simulation
% for each configuration...
for i_times=1:Nconfig
    
    disp('Simulation ' + string(i_times) + ' of ' + Nconfig)
    
    for j=1:Nparameters
        assignin('base',parameters.names{j}, parameters.values{j}(i_times));
    end
    
    % create initial conditions
    rng(seed,'twister'); % reproducible results
    x0Data=nan(Ntimes,N,2);
    for k_times=  1:Ntimes
        x0Data(k_times,:,:)=randCircle(N, initRadius);
    end
    
    tic
    parfor k_times=1:Ntimes
        x0=squeeze(x0Data(k_times,:,:));
        [T_r, success, final_e_theta, final_e_L, finalGRadial, finalGNormal, stopTime] = Simulator(x0, LinkNumber, G_radial, G_normal, regularity_thresh, compactness_thresh, Tmax, sigma, drawON, getMetrics, IntFunctionStruct, MaxSensingRadius, alpha, beta, dynamicLattice, AgentsRemoval);

        e_L(k_times)=final_e_L;
        e_theta(k_times)=final_e_theta;
        Tr_vec(k_times)=T_r;
        success_vec(k_times)=success;
        stopTime_vec(k_times) = stopTime;
        GNormal_vec(k_times) = finalGNormal;
    end
    
    %Average and Variance
    e_L_mean(i_times,:)=[mean(e_L), min(e_L), max(e_L)];
    e_theta_mean(i_times,:)=[mean(e_theta), min(e_theta), max(e_theta)];
    Tr_mean(i_times)=mean(Tr_vec,'omitnan');
    success_mean(i_times)=mean(success_vec);
    stopTime_max(i_times) = max(stopTime_vec,[], 'includenan');
    GNormal_mean(i_times,:)=[mean(GNormal_vec), min(GNormal_vec), max(GNormal_vec)];

    toc
end

cost_mean=[vecnorm([e_theta_mean(:,1)/regularity_thresh, e_L_mean(:,1)/compactness_thresh],2,2), vecnorm([e_theta_mean(:,2)/regularity_thresh, e_L_mean(:,2)/compactness_thresh],2,2), vecnorm([e_theta_mean(:,3)/regularity_thresh, e_L_mean(:,3)/compactness_thresh],2,2)];

%% Output in command window

fprintf('\n --- \nNtimes=%d seed=%d \n\nNconf \t|',Ntimes,seed)
for j_times=1:Nparameters
    fprintf('%s\t',string(parameters.names{j_times}));
end
fprintf('\t| e_t \t e_L \t T_r \t succ \t cost \n --- \n');

for i_times=1:Nconfig
    fprintf(' %d \t|', i_times)
    for j_times=1:Nparameters
        fprintf('%.2f\t\t',string(parameters.values{j_times}(i_times)));
    end
    fprintf('| %.2f \t %.2f \t %.2f \t %.2f \t %.2f \n',e_theta_mean(i_times),e_L_mean(i_times),Tr_mean(i_times),success_mean(i_times),cost_mean(i_times));
end
fprintf(' --- \nMean \t|')
for j_times=1:Nparameters
    fprintf('xxxx \t\t');
end
fprintf('| %.2f \t %.2f \t %.2f \t %.2f \t %.2f \n --- \n',mean(e_theta_mean(:,1)),mean(e_L_mean(:,1)),mean(Tr_mean),mean(success_mean),mean(cost_mean(:,1)));

%% Plots
% plot if Nparameters==1
if Nparameters==1
    
    if alpha==0 && beta==0
        figure %METRICS
        set(0, 'DefaultFigureRenderer', 'painters');
        subplot(2,1,1)
        title('Metrics','FontSize',16)
        hold on
        e_theta_line=plotWithShade(parameters.values{1}, e_theta_mean(:,1), e_theta_mean(:,2), e_theta_mean(:,3), 'b', 0.1);
        e_theta_star_line=yline(regularity_thresh, 'b--','LineWidth',2);
        axis([-inf inf 0 0.4])
        yticks([0.1 regularity_thresh 0.3 0.4 0.5 0.6])
        xticks(parameters.values{1})
        %legend([e_theta_line, e_theta_star_line],{'0','1'}, 'Interpreter','latex','FontSize',22)
        legend([e_theta_line, e_theta_star_line],{'$\bar{e}_{\theta}$','$e_{\theta}^*$'},'Interpreter','latex','FontSize',22)
        legend boxoff
        set(gca,'FontSize',14)
        box
        
        subplot(2,1,2)
        e_L_line=plotWithShade(parameters.values{1},e_L_mean(:,1),e_L_mean(:,2),e_L_mean(:,3), 'r', 0.1);
        e_L_star_line=yline(compactness_thresh,'r--','LineWidth',2);
        yticks([0.15 compactness_thresh 0.45 0.6 0.75 0.9])
        xticks(parameters.values{1})
        axis([-inf inf 0 0.6])
        %legend([e_L_line, e_L_star_line],{'2','3'},'Interpreter','latex','FontSize',22)
        legend([e_L_line, e_L_star_line],{'$\bar{e}_L$','$e_L^*$'},'Interpreter','latex','FontSize',22)
        legend boxoff
        xlabel(parameters.names{1})
        set(gca,'FontSize',14)
        box
    else
         figure %METRICS with Gain
         set(0, 'DefaultFigureRenderer', 'painters');
         subplot(3,1,1)
         title('Metrics','FontSize',16)
         hold on
         e_theta_line=plotWithShade(parameters.values{1}, e_theta_mean(:,1), e_theta_mean(:,2), e_theta_mean(:,3), 'b', 0.1);
         e_theta_star_line=yline(regularity_thresh, 'b--','LineWidth',2);
         axis([-inf inf 0 0.6])
         yticks([0 regularity_thresh 0.4 0.6])
         xticks(parameters.values{1})
         %legend([e_theta_line, e_theta_star_line],{'0','1'}, 'Interpreter','latex','FontSize',22)
         legend([e_theta_line, e_theta_star_line],{'$\bar{e}_{\theta}$','$e_{\theta}^*$'},'Interpreter','latex','FontSize',22)
         legend boxoff
         set(gca,'FontSize',14)
         box
         
         subplot(3,1,2)
         e_L_line=plotWithShade(parameters.values{1},e_L_mean(:,1),e_L_mean(:,2),e_L_mean(:,3), 'r', 0.1);
         e_L_star_line=yline(compactness_thresh,'r--','LineWidth',2);
         yticks([0 compactness_thresh 0.6 0.9])
         xticks(parameters.values{1})
         axis([-inf inf 0 0.6])
         %legend([e_L_line, e_L_star_line],{'2','3'},'Interpreter','latex','FontSize',22)
         legend([e_L_line, e_L_star_line],{'$\bar{e}_L$','$e_L^*$'},'Interpreter','latex','FontSize',22)
         legend boxoff
         set(gca,'FontSize',14)
         box
         
         subplot(3,1,3)
         line=plotWithShade(parameters.values{1},GNormal_mean(:,1),GNormal_mean(:,2),GNormal_mean(:,3), 'm', 0.1);
         xticks(parameters.values{1})
         yticks([0:2.5:7.5])
         axis([-inf inf 0 7.5])
         %legend([line],{'4'},'Interpreter','latex','FontSize',22)
         legend([line],{'$\bar{G}_n$'},'Interpreter','latex','FontSize',22)
         legend boxoff
         xlabel(parameters.names{1})
         set(gca,'FontSize',14)
         box
         set(gcf,'Position',[250 250 500 500])
    end
    
    figure %Success Rate and Cost
    subplot(2,1,1)
    set(0, 'DefaultFigureRenderer', 'painters');
    title('Success Rate','FontSize',16)
    hold on
    line=plot(parameters.values{1}, success_mean, '-ko','LineWidth',2);
    axis([-inf inf 0 1.1])
    legend([line],{'0'}, 'Interpreter','latex','FontSize',22)
    legend([line],{'Success rate'},'Interpreter','latex','FontSize',22)
    legend boxoff        
    set(gca,'FontSize',14)
    box
    
    subplot(2,1,2)
    line=plotWithShade(parameters.values{1},cost_mean(:,1),cost_mean(:,2),cost_mean(:,3), 'r', 0.1);
    level=yline(1,'r--','LineWidth',2);
    axis([-inf inf 0 2])
    legend([line, level],{'cost','1'},'Interpreter','latex','FontSize',22)
    legend boxoff
    xlabel(parameters.names{1})
    set(gca,'FontSize',14)
    box
end





