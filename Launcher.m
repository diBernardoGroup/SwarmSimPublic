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
outputDir='/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/stability of geometric lattices/simulations';


defaultParam;   % load default parameters

N=10;

avgSpeed0=1;
sigmaSpeed0=0.5;

Dynamics=struct('model','FirstOrder', 'sigma',0, 'vMax', 5);
%Dynamics=struct('model','SecondOrder', 'sigma',0.1, 'vMax', inf);
%Dynamics=struct('model','CoupledSDEs', 'rateSpeed', 1, 'avgSpeed', avgSpeed0, 'sigmaSpeed', 1, 'rateOmega', 1, 'sigmaOmega', @(x)2*max(1-x/3,0), 'omega', zeros(N,1));
%Dynamics=struct('model','LevyWalk', 'alpha',0.005, 'sigma', 0.25);

smoothing = false;      

%% Create Initial Conditions
%rng(1,'twister'); % set the randomn seed to have reproducible results

%x0=randCircle(N, 2);                 % initial conditions drawn from a uniform disc
%x0 = normrnd(0,0.1*sqrt(N),N,2);    % initial conditions drawn from a normal distribution
%x0 = perfectLactice(N, LinkNumber); % initial conditions on a correct lattice
x0 = perfectLactice(N, LinkNumber) + randCircle(N, 0.6); % initial conditions on a deformed lattice

speeds0 = abs(normrnd(avgSpeed0,sigmaSpeed0,N,1));
theta0 = 2*pi*rand(N,1)-pi;
v0 = speeds0 .* [cos(theta0), sin(theta0)];
%v0 = zeros(size(x0));

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

%% PLOTS
close all

% create folder, save data and parameters
counter=1;
while exist(fullfile(outputDir,[Dynamics.model,num2str(counter)]),'dir')
    counter=counter+1;
end
path=fullfile(outputDir, [Dynamics.model,num2str(counter)]);
mkdir(path)
save(fullfile(path, 'data'))

fileID = fopen(fullfile(path, 'parameters.txt'),'wt');
fprintf(fileID,'Date: %s\n',datestr(now, 'dd/mm/yy'));
fprintf(fileID,'Time: %s\n\n',datestr(now, 'HH:MM'));
fprintf(fileID,'Parameters:\n\n');
fprintf(fileID,'N= %d\n',N);
fprintStruct(fileID,Simulation)
fprintStruct(fileID,Dynamics)
fprintf(fileID,'smoothing= %s\n',mat2str(smoothing));
fclose(fileID);

% TRAJECTORIES
figure
if Simulation.drawTraj; plotTrajectory(xVec, false, [0,0.7,0.9]); end
if isfield(LocalIntFunction, 'DistanceRange')
    plotSwarmInit(squeeze(xVec(end,:,:)), Simulation.Tmax, LocalIntFunction.DistanceRange(1), LocalIntFunction.DistanceRange(2));
else
    plotSwarmInit(squeeze(xVec(end,:,:)), Simulation.Tmax, inf, inf);
end
saveas(gcf, fullfile(path, 'trajectories'))
saveas(gcf, fullfile(path, 'trajectories'),'png')

if ~strcmp(GlobalIntFunction.function,'None') % RADIAL INTERACTION FUNCTION
    figure 
    hold on
    fplot(@(x) RadialInteractionForce(x, GlobalIntFunction),[0, 3])
    plot([1], [0], 'r.','MarkerSize', 40)
    yticks([-1 0 1])
    xticks([0:3])
    ylim([-0.2 1.2])
    grid on
    title('f_r(d)')
    xlabel('d')
    set(gca,'FontSize',14)
    saveas(gcf,fullfile(path, 'radial_inter_func'))
    saveas(gcf,fullfile(path, 'radial_inter_func'),'png')
end

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
saveas(gcf,fullfile(path, 'time_plot'))
saveas(gcf,fullfile(path, 'time_plot'),'png')

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
% saveas(gcf,fullfile(path, 'scatter_plot'))
% saveas(gcf,fullfile(path, 'scatter_plot'),'png')
% 
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
% saveas(gcf,fullfile(path, 'scatter_plot_normalised'))
% saveas(gcf,fullfile(path, 'scatter_plot_normalised'),'png')
