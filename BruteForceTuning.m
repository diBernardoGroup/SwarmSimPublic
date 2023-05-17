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

Ntimes=1;       % How many times simulator is launched on each config (couple of gains)

D=2;            %number of dimensions [2 or 3]

LinkNumber=4;   % number of links (6=triangular lattice, 4=square lattice, 3=hexagonal lattice) (L)

defaultParam;   % load default parameters

window=1;       % smoothing window size (1 = no smoothing)

%ranges of values of control gains
G_r_min = 0;
G_r_max = 10;
G_n_min = 0;
G_n_max = 10;

%grid of values
step=10;
G_r_vec = [G_r_min:step:G_r_max];
G_n_vec = [G_n_min:step:G_n_max];

%% Preallocate variables

m       = nan(length(G_r_vec), length(G_n_vec), Ntimes);
rigid   = nan(length(G_r_vec), length(G_n_vec), Ntimes);
e_d_max = nan(length(G_r_vec), length(G_n_vec), Ntimes);

% e_theta = zeros(1, Ntimes);
% e_L = zeros(1, Ntimes);
% e_d = zeros(1, Ntimes);
% Tr_vec = zeros(1, Ntimes);
% success_vec = zeros(1, Ntimes);
% stopTime_vec = zeros(1, Ntimes);
% e_theta_mean = nan(length(G_r_vec), length(G_n_vec), Ntimes);
% e_L_mean = zeros(length(G_r_vec), length(G_n_vec));
% e_d_mean = zeros(length(G_r_vec), length(G_n_vec));
% Tr_mean = zeros(length(G_r_vec), length(G_n_vec));
% success_mean = zeros(length(G_r_vec), length(G_n_vec));
% stopTime_mean = zeros(length(G_r_vec), length(G_n_vec));
% stopTime_max = zeros(length(G_r_vec), length(G_n_vec));

%% Create initial conditions
rng(0,'twister'); % reproducible results
x0Data=nan(Ntimes,N,2);
for k_times=1:Ntimes
    x0Data(k_times,:,:)=randCircle(N, 2, D); % circle of radius 2
end

%% Simulate Ntimes over all values of control gains
for i_times= 1:length(G_r_vec)
    for j_times= 1:length(G_n_vec)
        GlobalIntFunction.Gain  = G_r_vec(i_times); 
        LocalIntFunction.Gain   = G_n_vec(j_times);
        disp('Simulations batch ' + string((i_times-1)*length(G_n_vec)+j_times) + ' of ' + string(length(G_r_vec)*length(G_n_vec)))

        tic
        for k_times=  1:Ntimes
            x0=squeeze(x0Data(k_times,:,:));
            v0 = zeros(N,D);
                    
            %Simulation
            [xVec, uVec] = Simulator(x0, v0, Simulation, Dynamics, GlobalIntFunction, LocalIntFunction);
            
            %Analysis of the terminal configuration
            x_f=squeeze(xVec(end,:,:));
            B = buildIncidenceMatrix(x_f, Rmax);
            m(i_times,j_times,k_times)=size(B,2);
            M = buildRigidityMatrix(x_f, B);
            rigid(i_times,j_times,k_times) = rank(M)==D*N-D*(D+1)/2;
            e_d_max(i_times,j_times,k_times) = getMaxLinkLengthError(x_f, 1, 0, Rmax);
        end
        toc
        
%         e_theta_mean(i_times, j_times) = mean(e_theta);
%         e_L_mean(i_times, j_times) = mean(e_L);
%         e_d_mean(i_times, j_times) = mean(e_d);
%         Tr_mean(i_times, j_times) = mean(Tr_vec, 'omitnan');
%         success_mean(i_times, j_times) = mean(success_vec);
%         stopTime_mean(i_times, j_times) = mean(stopTime_vec, 'omitnan');
%         stopTime_max(i_times, j_times) = max(stopTime_vec,[], 'includenan');
    end
end

% average over the initial conditions
m_mean = mean(e_d_max,3);
rigid_mean = mean(e_d_max,3);
e_d_max_mean = mean(e_d_max,3);

% smoothing over the parameters space
m_mean=movmean2(m_mean,window);
rigid_mean=movmean2(rigid_mean,window);
e_d_max_mean=movmean2(e_d_max_mean,window);


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


