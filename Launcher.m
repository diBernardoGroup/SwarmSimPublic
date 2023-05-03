%
%Launcher Set the parameters and launch a single simulation of the swarm.
%   Also robustness tests can be run, see AgentsRemoval, NoiseTest and dynamicLattice
%
%   See also: BruteForceTuning, SequentialLauncher, StabilityAnalysis
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%

%% Clear environment
close all
clear

%% Parameters

defaultParam;   % load default parameters

D=3; %number of dimensions [2 or 3]
N=20;

% avgSpeed0=1;
% sigmaSpeed0=0.5;

Dynamics=struct('model','FirstOrder', 'sigma',0, 'vMax', inf);
%Dynamics=struct('model','SecondOrder', 'sigma',0.1, 'vMax', inf);
%Dynamics=struct('model','CoupledSDEs', 'rateSpeed', 1, 'avgSpeed', avgSpeed0, 'sigmaSpeed', 1, 'rateOmega', 1, 'sigmaOmega', @(x)2*max(1-x/3,0), 'omega', zeros(N,1));
%Dynamics=struct('model','LevyWalk', 'alpha',0.005, 'sigma', 0.25);

smoothing = false;

%% Create Initial Conditions
%rng(1,'twister'); % set the randomn seed to have reproducible results

%x0=randCircle(N, 2, D);                 % initial conditions drawn from a uniform disc
%x0 = normrnd(0,0.1*sqrt(N),N,D);    % initial conditions drawn from a normal distribution
x0 = perfectLactice(N, LinkNumber, D); % initial conditions on a correct lattice
%x0 = perfectLactice(N, LinkNumber, D) + randCircle(N, delta, D); % initial conditions on a deformed lattice

% speeds0 = abs(normrnd(avgSpeed0,sigmaSpeed0,N,1));
% theta0 = 2*pi*rand(N,1)-pi;
% v0 = speeds0 .* [cos(theta0), sin(theta0)];
v0 = zeros(size(x0));

%% Run Simulation
[xVec] = Simulator(x0, v0, Simulation, Dynamics, GlobalIntFunction, LocalIntFunction);

%% Analysis
if smoothing
    xVec = movmean(xVec,3);
    xVec = movmean(xVec,3);
end

timeInstants = 0:Simulation.deltaT:Simulation.Tmax;

% derivate quantities
[~, vVec] = gradient(xVec, 1, Simulation.deltaT, 1);
speed = vecnorm(vVec,2,3);
theta = atan2(vVec(:,:,2), vVec(:,:,1));
for i=1:length(timeInstants)-1
    omega(i,:) = angleBetweenVectors(squeeze(vVec(i,:,:)),squeeze(vVec(i+1,:,:)))';
end
omega(length(timeInstants),:) = angleBetweenVectors(squeeze(vVec(length(timeInstants)-1,:,:)),squeeze(vVec(length(timeInstants),:,:)))';
omega=omega/Simulation.deltaT;

% metrics
for i=1:length(timeInstants) % for each time instant...
    x=squeeze(xVec(i,:,:));
    
    e_d(i) = getAvgLinkLengthError(x, 1, 0, Rmax);          % avg distance from the deisred link length
    e_d_max(i) = getMaxLinkLengthError(x, 1, 0, Rmax);      % max distance from the deisred link length.
    % e_d_max<=(Rmax-1)preserves all the links.
    % e_d_max(0)<= 2*delta
    
    B = buildIncidenceMatrix(x, Rmax);                      % incidence matrix
    m(i)=size(B,2);                                         % number of links
    M = buildRigidityMatrix(x, B);                          % rigidity matrix
    
    rigidity(i) = rank(M)==D*N-D*(D+1)/2;                   % check infinitesimal rigidity
end

%% PLOTS
close all

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
    fprintf(fileID,'Launcher\n\n');
    fprintf(fileID,'Date: %s\n',datestr(now, 'dd/mm/yy'));
    fprintf(fileID,'Time: %s\n\n',datestr(now, 'HH:MM'));
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
figure
if isfield(LocalIntFunction, 'DistanceRange')
    plotSwarmInit(squeeze(xVec(end,:,:)), Simulation.Tmax, LocalIntFunction.DistanceRange(1), LocalIntFunction.DistanceRange(2));
else
    plotSwarmInit(squeeze(xVec(end,:,:)), Simulation.Tmax, inf, inf);
end
if Simulation.drawTraj; plotTrajectory(xVec, false, [0,0.7,0.9]); end
if outputDir
    saveas(gcf, fullfile(path, 'trajectories'))
    saveas(gcf, fullfile(path, 'trajectories'),'png')
end

% if ~strcmp(GlobalIntFunction.function,'None') % RADIAL INTERACTION FUNCTION
%     figure
%     hold on
%     fplot(@(x) RadialInteractionForce(x, GlobalIntFunction),[0, 3])
%     plot([1], [0], 'r.','MarkerSize', 40)
%     yticks([-1 0 1])
%     xticks([0:3])
%     ylim([-0.2 1.2])
%     grid on
%     title('f_r(d)')
%     xlabel('d')
%     set(gca,'FontSize',14)
%     if outputDir
%     saveas(gcf,fullfile(path, 'radial_inter_func'))
%     saveas(gcf,fullfile(path, 'radial_inter_func'),'png')
%     end
% end

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


figure % SPEED and ANGULAR VELOCITY
subplot(2,4,[1 2 3])
plotWithShade(timeInstants, median(speed,2), min(speed, [], 2), max(speed, [], 2), 'b', 0.3);
xlabel('t [s]')
ylabel('speed')
rng=ylim;
box on
subplot(2,4,4)
h=histogram(speed(:),'Orientation','horizontal');
ylim(rng);
set(gca,'xtick',[])
subplot(2,4,[5 6 7])
plotWithShade(timeInstants, median(abs(omega),2), min(abs(omega), [], 2), max(abs(omega), [], 2), 'b', 0.3);
xlabel('t [s]')
ylabel('ang. vel. [rad/s]')
rng=ylim;
box on
subplot(2,4,8)
h=histogram(abs(omega(:)),'Orientation','horizontal');
ylim(rng);
set(gca,'xtick',[])
if outputDir
    saveas(gcf,fullfile(path, 'time_plot'))
    saveas(gcf,fullfile(path, 'time_plot'),'png')
end

% figure % SCATTER PLOT - SPEED and ANGULAR VELOCITY
% s=scatterhist(speed(:),abs(omega(:)), 'Location','NorthEast','Direction','out');
% xlabel(s,'speed')
% ylabel(s,'ang. vel. [rad/s]')
% s(1).YAxisLocation = 'left';
% s(1).XAxisLocation = 'bottom';
% s(2).Position = [0.1    0.82   0.7    0.125];
% s(3).Position = [0.82   0.1    0.125    0.7];
% s(1).Position(3) = 0.7;
% s(1).Position(4) = 0.7;
% if outputDir
% saveas(gcf,fullfile(path, 'scatter_plot'))
% saveas(gcf,fullfile(path, 'scatter_plot'),'png')
% end

% figure % SCATTER PLOT - NORMALIZED SPEED and ANGULAR VELOCITY
% normalized_speed=speed./median(speed);
% s=scatterhist(normalized_speed(:),abs(omega(:)), 'Location','NorthEast','Direction','out');
% xlabel(s,'speed / agent median speed')
% ylabel(s,'ang. vel. [rad/s]')
% s(1).YAxisLocation = 'left';
% s(1).XAxisLocation = 'bottom';
% s(2).Position = [0.1    0.82   0.7    0.125];
% s(3).Position = [0.82   0.1    0.125    0.7];
% s(1).Position(3) = 0.7;
% s(1).Position(4) = 0.7;
% if outputDir
% saveas(gcf,fullfile(path, 'scatter_plot_normalised'))
% saveas(gcf,fullfile(path, 'scatter_plot_normalised'),'png')
% end

figure % e_d_max
set(gca,'FontSize',14)
set(0, 'DefaultFigureRenderer', 'painters');
set(gcf,'Position',[100 100 560 420*0.6])
hold on
line=plot(timeInstants, e_d_max, 'b');
yline(Rmax-1,'--','LineWidth',2)
yticks([0:0.1:0.3, Rmax-1, 0.4:0.1:1])
set(gca,'YTickLabel',{[0:0.1:0.3], 'R_a-R', [0.4:0.1:1]})
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
%     plot(timeInstants,m)
%     title('m', 'Interpreter','latex','FontSize',22)
%     xlabel('t', 'Interpreter','latex','FontSize',22)
%     set(gca,'FontSize',14)
%     box
%     grid

figure % rigidity
set(gca,'FontSize',14)
set(gcf,'Position',[100 100 560 420*0.6])
hold on
plot(timeInstants,rigidity)
axis([-inf inf -0.05 1.05])
title('$\rho$', 'Interpreter','latex','FontSize',22)
xlabel('t', 'Interpreter','latex','FontSize',22)
box
grid
if outputDir
    saveas(gcf,fullfile(path, 'rigidity'))
    saveas(gcf,fullfile(path, 'rigidity'),'png')
end