%
%StabilityAnalysis Set the parameters and launch multiple simulations from different initial conditions.
%   Also robustness tests can be run, see AgentsRemoval, NoiseTest and dynamicLattice
%
%   See also: Launcher, SequentialLauncher
%
%   Authors:    Andrea Giusti
%   Date:       2022
%

%% Clear environment
close all
clear
clc

%% Parameters

Ntimes=3;              % How many simulations are launched for each configuration

defaultParam;           % load default parameters

D=3; %number of dimensions [2 or 3]
N=20;

% avgSpeed0=1;
% sigmaSpeed0=0.5;

Dynamics=struct('model','FirstOrder', 'sigma',0, 'vMax', inf);
%Dynamics=struct('model','SecondOrder', 'sigma',0.1, 'vMax', inf);
%Dynamics=struct('model','CoupledSDEs', 'rateSpeed', 1, 'avgSpeed', avgSpeed0, 'sigmaSpeed', 1, 'rateOmega', 1, 'sigmaOmega', @(x)2*max(1-x/3,0), 'omega', zeros(N,1));
%Dynamics=struct('model','LevyWalk', 'alpha',0.005, 'sigma', 0.25);

smoothing = false;

rng(0,'twister');       % set the randomn seed to have reproducible results


%% Run Simulation
disp(['Running ',num2str(Ntimes),' simulations:'])
for rep=1:Ntimes
    %% Create Initial Conditions
    %x0=randCircle(N, 2, D);                               % initial conditions drawn from a uniform disc
    %x0 = normrnd(0,0.1*sqrt(N),N,2);                   % initial conditions drawn from a normal distribution
    %x0 = perfectLactice(N, LinkNumber, true, true, (sqrt(N)+1)^2);        % initial conditions on a correct lattice
    x0 = perfectLactice(N, LinkNumber, D) + randCircle(N, delta, D); % initial conditions on a deformed lattice
    
    v0 = zeros(size(x0));
    
    %% Run Simulation
    [xVec(rep,:,:,:)] = Simulator(x0, v0, Simulation, Dynamics, GlobalIntFunction, LocalIntFunction);
    
    %% ANALYSIS
    
    timeInstants = 0:Simulation.deltaT:Simulation.Tmax;
    
    % metrics
    for i=1:length(timeInstants) % for each time instant...
        x=squeeze(xVec(rep,i,:,:));
        
        e_d(rep,i) = getAvgLinkLengthError(x, 1, 0, Rmax);          % avg distance from the deisred link length
        e_d_max(rep,i) = getMaxLinkLengthError(x, 1, 0, Rmax);      % max distance from the deisred link length.
        % e_d_max<=(Rmax-1)preserves all the links.
        % e_d_max(0)<= 2*delta
        
        B = buildIncidenceMatrix(x, Rmax);                      % incidence matrix
        m(rep,i)=size(B,2);                                         % number of links
        M = buildRigidityMatrix(x, B);                          % rigidity matrix
        
        rigidity(rep,i) = rank(M)==D*N-D*(D+1)/2;                   % check infinitesimal rigidity
    end
end

%% PLOTS

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
    fprintf(fileID,'StabilityAnalysis\n\n');
    fprintf(fileID,'Date: %s\n',datestr(now, 'dd/mm/yy'));
    fprintf(fileID,'Time: %s\n\n',datestr(now, 'HH:MM'));
    fprintf(fileID,'Ntimes= %d\n\n',Ntimes);
    fprintf(fileID,'Parameters:\n\n');
    fprintf(fileID,'N= %d\n',N);
    fprintf(fileID,'D= %d\n',size(x0,2));
    fprintStruct(fileID,Simulation)
    fprintf(fileID,'Dynamics:\n');
    fprintStruct(fileID,Dynamics)
    fprintf(fileID,'GlobalIntFunction:\n');
    fprintStruct(fileID,GlobalIntFunction)
    fprintf(fileID,'LocalIntFunction:\n');
    fprintStruct(fileID,LocalIntFunction)
    fprintf(fileID,'smoothing= %s\n',mat2str(smoothing));
    fclose(fileID);
end

% SWARM
if isfield(LocalIntFunction, 'DistanceRange')
    figure
    plotSwarmInit(squeeze(xVec(1,1,:,:)), 0, LocalIntFunction.DistanceRange(1), LocalIntFunction.DistanceRange(2))
    figure
    plotSwarmInit(squeeze(xVec(1,ceil(length(timeInstants)/2),:,:)), Simulation.Tmax/2, LocalIntFunction.DistanceRange(1), LocalIntFunction.DistanceRange(2))
    figure
    plotSwarmInit(squeeze(xVec(1,length(timeInstants),:,:)), Simulation.Tmax, LocalIntFunction.DistanceRange(1), LocalIntFunction.DistanceRange(2))
else
    figure
    plotSwarmInit(squeeze(xVec(1,1,:,:)), 0, inf, inf)
    figure
    plotSwarmInit(squeeze(xVec(1,ceil(length(timeInstants)/2),:,:)), Simulation.Tmax/2, inf, inf)
    figure
    plotSwarmInit(squeeze(xVec(1,length(timeInstants),:,:)), Simulation.Tmax, inf, inf)
end

%     figure % RADIAL INTERACTION FUNCTION
%     hold on
%     set(gca,'FontSize',14)
%     fplot(@(x) RadialInteractionForce(x, RadialIntFunction),[0, 2], 'LineWidth', 2)
%     plot([1], [0], 'r.','MarkerSize', 30)
%     yticks([-1:0.5:1])
%     xticks([0:0.5:1, Rmax, 1.5:0.5:3])
%     set(gca,'XTickLabel',{[0:0.5:1], 'R_a', [1.5:0.5:3]})
%     ylim([-0.2 1.2])
%     grid on
%     title('$f(z)$','Interpreter','latex', 'FontSize', 22)
%     xlabel('$z$','Interpreter','latex', 'FontSize', 22)
%     box on
%     grid on

%     figure % NORMAL INTERACTION FORCE
%     hold on
%     fplot(@(alfa) NormalInteractionForce(alfa, LinkNumber),[-pi/LinkNumber, pi/LinkNumber])
%     plot([0], [0], 'r.','MarkerSize', 30)
%     ylim([-1.2 1.2])
%     xlim([-pi/LinkNumber pi/LinkNumber])
%     yticks([-1 0 1])
%     xticks([-pi/LinkNumber, 0, pi/LinkNumber])
%     set(gca,'XTickLabel',{'-\pi/4','0','\pi/4'})
%     grid on
%     title('f_n(\theta)')
%     xlabel('\theta')
%     set(gca,'FontSize',14)

figure % e_d_max
set(gca,'FontSize',14)
set(0, 'DefaultFigureRenderer', 'painters');
set(gcf,'Position',[100 100 560 420*0.6])
hold on
line=plotWithShade(timeInstants, mean(e_d_max), min(e_d_max), max(e_d_max), 'b', 0.2);
yline(Rmax-1,'--','LineWidth',2)
yticks(sort([0:0.1:1, Rmax-1]))
%set(gca,'YTickLabel',{[0:0.1:0.3], 'R_a-R', [0.4:0.1:1]})
%title('$e_{d,max}$', 'Interpreter','latex','FontSize',22)
%title('$\max_{i\in\mathcal{E}} |\Vert \mathbf{r}_{i} \Vert - R |$', 'Interpreter','latex','FontSize',22)
%legend([line],{'$e$'},'Interpreter','latex','FontSize',22)
ylabel('$e$', 'Interpreter','latex','FontSize',22, 'rotation',0,'VerticalAlignment','middle')
xlabel('t', 'Interpreter','latex','FontSize',22)
box
grid
if outputDir
    saveas(gcf,fullfile(path, 'e_d_max'))
    saveas(gcf,fullfile(path, 'e_d_max'),'png')
end


%     figure % m
%     plot(linspace(0,min(stopTime,Tmax,'omitnan'),time_instants),m)
%     title('m', 'Interpreter','latex','FontSize',22)
%     xlabel('t', 'Interpreter','latex','FontSize',22)
%     set(gca,'FontSize',14)
%     box
%     grid

figure % rigidity
set(gca,'FontSize',14)
set(gcf,'Position',[100 100 560 420*0.6])
hold on
plot(timeInstants,mean(rigidity),'r')
axis([-inf inf -0.05 1.05])
title('$\rho$', 'Interpreter','latex','FontSize',22)
xlabel('t', 'Interpreter','latex','FontSize',22)
box
grid
if outputDir
    saveas(gcf,fullfile(path, 'rigidity'))
    saveas(gcf,fullfile(path, 'rigidity'),'png')
end