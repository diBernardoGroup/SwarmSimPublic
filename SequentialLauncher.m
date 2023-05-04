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

Ntimes=3;              % How many simulations are launched for each configuration

defaultParam;           % load default parameters

D=2; %number of dimensions [2 or 3]
N=20;

seed=0;        % seed for random generator, if negative it is not set

smoothing = false;

%% variable parameters
% one or multiple parameters can be modified at the same time
% all the modified parameters must have the same number of values
% the values specified here overwrite the default ones
% IntFunctionStruct cannot be modified
% parameters.names={'initRadius'};
% parameters.values={[sqrt(25/25) sqrt(50/25) sqrt(100/25)]};
parameters.names={'delta'};
parameters.values={[0:0.25:1]};

%% Preallocate
Nparameters=length(parameters.names);
Nconfig=length(parameters.values{1});

e_L=nan(Nconfig,Ntimes);
e_theta=nan(Nconfig,Ntimes);
Tr_vec=nan(Nconfig,Ntimes);
success_vec=nan(Nconfig,Ntimes);
stopTime_vec = nan(Nconfig, Ntimes);
GNormal_vec = nan(Nconfig, Ntimes);
e_d_max_vec = nan(Nconfig, Ntimes);
V_vec = nan(Nconfig, Ntimes);
rigid_vec = nan(Nconfig, Ntimes);

%% Simulation
% for each configuration...
for i_times=1:Nconfig
    
    disp('Simulation ' + string(i_times) + ' of ' + Nconfig + ':')
    
    for j=1:Nparameters
        assignin('base',parameters.names{j}, parameters.values{j}(i_times));
    end
    
    % create initial conditions
    if seed>=0
        rng(seed,'twister'); % reproducible results
    end
    x0Data=nan(Ntimes,N,2);
    for k_times=  1:Ntimes
        %x0Data(k_times,:,:)=randCircle(N, 2, D);                               % initial conditions drawn from a uniform disc
        %x0Data(k_times,:,:) = normrnd(0,0.1*sqrt(N),N,2);                   % initial conditions drawn from a normal distribution
        %x0Data(k_times,:,:) = perfectLactice(N, LinkNumber, true, true, (sqrt(N)+1)^2);        % initial conditions on a correct lattice
        x0Data(k_times,:,:) = perfectLactice(N, LinkNumber, D) + randCircle(N, delta, D); % initial conditions on a deformed lattice
        
        v0 = zeros(N,D);
    end
    
    tic
    for k_times=1:Ntimes
        x0=squeeze(x0Data(k_times,:,:));
        
        [xVec(k_times,:,:,:)] = Simulator(x0, v0, Simulation, Dynamics, GlobalIntFunction, LocalIntFunction);
        
        %         e_L(k_times)=final_e_L;
        %         e_theta(k_times)=final_e_theta;
        %         if(size(T_r)==[1,2]); Tr_vec(k_times)=T_r; end
        %         success_vec(k_times)=success;
        %         stopTime_vec(k_times) = stopTime;
        %         GNormal_vec(k_times) = finalGNormal;
        
        x=squeeze(xVec(k_times,end,:,:));
        
        B = buildIncidenceMatrix(x, Rmax);
        m(i_times,k_times)=size(B,2);
        M = buildRigidityMatrix(x, B);
        rigid_vec(i_times,k_times) = rank(M)==2*N-3;
        e_d_max_vec(i_times,k_times) = getMaxLinkLengthError(x, 1, 0, Rmax);   % max distance from the deisred link length.
    end
    toc
end

cost=[vecnorm([mean(e_theta,2)/regularity_thresh, mean(e_L,2)/compactness_thresh],2,2), vecnorm([ min(e_theta,[],2)/regularity_thresh,  min(e_L,[],2)/compactness_thresh],2,2), vecnorm([ max(e_theta,[],2)/regularity_thresh,  max(e_L,[],2)/compactness_thresh],2,2)];

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
    fprintf('| %.2f \t %.2f \t %.2f \t %.2f \t %.2f \n',mean(e_theta(i_times,:)),mean(e_L(i_times,:)),mean(Tr_vec(i_times,:)),mean(success_vec(i_times,:)),mean(cost(i_times,:)));
end
fprintf(' --- \nMean \t|')
for j_times=1:Nparameters
    fprintf('xxxx \t\t');
end
fprintf('| %.2f \t %.2f \t %.2f \t %.2f \t %.2f \n --- \n',mean(e_theta,'all'),mean(e_L,'all'),mean(Tr_vec,'all'),mean(success_vec,'all'),mean(cost(:,1)));

%% Plots

% create folder, save data and parameters
if outputDir
    counter=1;
    while exist(fullfile(outputDir,[datestr(now, 'yyyy_mm_dd_'),Dynamics.model,'_',num2str(counter)]),'dir')
        counter=counter+1;
    end
    path=fullfile(outputDir, [datestr(now, 'yyyy_mm_dd_'),Dynamics.model,'_',num2str(counter)]);
    mkdir(path)
    save(fullfile(path, 'data'))
    
    fileID = fopen(fullfile(path, 'parameters.txt'),'wt');
    fprintf(fileID,'SequentialLauncher\n\n');
    fprintf(fileID,'Date: %s\n',datestr(now, 'dd/mm/yy'));
    fprintf(fileID,'Time: %s\n\n',datestr(now, 'HH:MM'));
    fprintf(fileID,'Ntimes= %d\n\n',Ntimes);
    fprintf(fileID,'Parameters:\n\n');
    fprintf(fileID,'N= %d\n',N);
    fprintf(fileID,'D= %d\n\n',size(x0,2));
    fprintf(fileID,'Simulation parameters:\n');
    fprintStruct(fileID,Simulation)
    fprintf(fileID,'Changing parameters:\n');
    fprintStruct(fileID,parameters)
    fprintf(fileID,'Dynamics:\n');
    fprintStruct(fileID,Dynamics)
    fprintf(fileID,'GlobalIntFunction:\n');
    fprintStruct(fileID,GlobalIntFunction)
    fprintf(fileID,'LocalIntFunction:\n');
    fprintStruct(fileID,LocalIntFunction)
    fprintf(fileID,'smoothing= %s\n',mat2str(smoothing));
    fclose(fileID);
end

% plot if Nparameters==1
if Nparameters==1
    
    if alpha==0 && beta==0
        %         figure %METRICS
        %         set(0, 'DefaultFigureRenderer', 'painters');
        %         subplot(2,1,1)
        %         title('Metrics','FontSize',16)
        %         hold on
        %         %e_theta_line=plotWithShade(parameters.values{1}, e_theta_mean(:,1), e_theta_mean(:,2), e_theta_mean(:,3), 'b', 0.1);
        %         e_theta_line=plotWithShade(parameters.values{1},  mean(e_theta,2), min(e_theta,[],2), max(e_theta,[],2), 'b', 0.1);
        %         e_theta_star_line=yline(regularity_thresh, 'b--','LineWidth',2);
        %         axis([-inf inf 0 0.4])
        %         yticks([0.1 regularity_thresh 0.3 0.4 0.5 0.6])
        %         xticks(parameters.values{1})
        %         %legend([e_theta_line, e_theta_star_line],{'0','1'}, 'Interpreter','latex','FontSize',22)
        %         legend([e_theta_line, e_theta_star_line],{'$\bar{e}_{\theta}$','$e_{\theta}^*$'},'Interpreter','latex','FontSize',22)
        %         legend boxoff
        %         set(gca,'FontSize',14)
        %         box
        %
        %         subplot(2,1,2)
        %         %e_L_line=plotWithShade(parameters.values{1},e_L_mean(:,1),e_L_mean(:,2),e_L_mean(:,3), 'r', 0.1);
        %         e_L_line=plotWithShade(parameters.values{1}, mean(e_L,2), min(e_L,[],2), max(e_L,[],2), 'r', 0.1);
        %         e_L_star_line=yline(compactness_thresh,'r--','LineWidth',2);
        %         yticks([0.15 compactness_thresh 0.45 0.6 0.75 0.9])
        %         xticks(parameters.values{1})
        %         axis([-inf inf 0 0.6])
        %         %legend([e_L_line, e_L_star_line],{'2','3'},'Interpreter','latex','FontSize',22)
        %         legend([e_L_line, e_L_star_line],{'$\bar{e}_L$','$e_L^*$'},'Interpreter','latex','FontSize',22)
        %         legend boxoff
        %         xlabel(parameters.names{1})
        %         set(gca,'FontSize',14)
        %         box
        
        figure %e_d_max and rigidity
        set(0, 'DefaultFigureRenderer', 'painters');
        subplot(2,1,1)
        hold on
        line=plotWithShade(parameters.values{1}, mean(e_d_max_vec,2), min(e_d_max_vec,[],2), max(e_d_max_vec,[],2), 'b', 0.1); %e_d_max_mean(:,1),e_d_max_mean(:,2),e_d_max_mean(:,3), 'b', 0.1);
        xticks(parameters.values{1})
        set(gca,'FontSize',14)
        ylabel('$e$', 'Interpreter','latex','FontSize',22, 'rotation',0,'VerticalAlignment','middle')
        %legend([line],{'$e$'},'Interpreter','latex','FontSize',22)
        %legend boxoff
        box on
        grid
        
        subplot(2,1,2)
        rigidity_line=plot(parameters.values{1}, mean(rigid_vec,2),'Marker','o','Color','r','LineWidth',2,'MarkerSize',6);
        xticks(parameters.values{1})
        xlabel(parameters.names{1})
        set(gca,'FontSize',14)
        ylabel('$\rho$', 'Interpreter','latex','FontSize',22, 'rotation',0,'VerticalAlignment','middle')
        %legend([rigidity_line],{'$\rho$'},'Interpreter','latex','FontSize',22)
        %legend boxoff
        box on
        grid
        if outputDir
            saveas(gcf,fullfile(path, 'e_rho'))
            saveas(gcf,fullfile(path, 'e_rho'),'png')
        end
        
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
    
    %     figure %Success Rate and Cost
    %     subplot(2,1,1)
    %     set(0, 'DefaultFigureRenderer', 'painters');
    %     title('Success Rate','FontSize',16)
    %     hold on
    %     line=plot(parameters.values{1}, mean(success_vec,2), '-ko','LineWidth',2);
    %     axis([-inf inf 0 1.1])
    %     legend([line],{'0'}, 'Interpreter','latex','FontSize',22)
    %     legend([line],{'Success rate'},'Interpreter','latex','FontSize',22)
    %     legend boxoff
    %     set(gca,'FontSize',14)
    %     box
    %
    %     subplot(2,1,2)
    %     line=plotWithShade(parameters.values{1},cost(:,1),cost(:,2),cost(:,3), 'r', 0.1);
    %     level=yline(1,'r--','LineWidth',2);
    %     axis([-inf inf 0 2])
    %     legend([line, level],{'cost','1'},'Interpreter','latex','FontSize',22)
    %     legend boxoff
    %     xlabel(parameters.names{1})
    %     set(gca,'FontSize',14)
    %     box
end





