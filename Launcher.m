%
%Launcher Set the parameters and launch a single simulation of the swarm.
%   Also robustness tests can be run, see AgentsRemoval, NoiseTest and dynamicLattice
%
%   See also: MultiLauncher, SequentialLauncher
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%

%% Clear environment
close all
clear

%% Parameters

D=2;                        % number of dimensions [2 or 3]

defaultParam;               % load default parameters

identification=readtable('/Volumes/DOMEPEN/Experiments/2023_06_15_Euglena_7/tracking_2023_10_16/identification.txt');
ids=randsample(length(identification.agents),N, true, ones(length(identification.agents),1));
agents = identification(ids,:);
Dynamics=struct('model','IndependentSDEsWithInput', ...
    'avgSpeed',agents.mu_s, 'rateSpeed', agents.theta_s, 'sigmaSpeed', agents.sigma_s, 'gainSpeed', agents.alpha_s,...
    'rateOmega', agents.theta_w, 'sigmaOmega', agents.sigma_w, 'gainOmega', agents.alpha_w, 'omega', normrnd(0,agents.std_w,N,1));
Environment = struct();
Environment.EnvUniform.Times  = [0 5 5.5 25 25.5 40]; 
Environment.EnvUniform.Values = [0 0   1  1    0  0]; 

Simulation.drawON=true;     % draw swarm during simulation (if N is large slows down the simulation)

%% Create Initial Conditions
%rng(1,'twister'); % set the randomn seed to have reproducible results

x0=randCircle(N, 1000, D);                 % initial conditions drawn from a uniform disc
%x0 = normrnd(0,0.1*sqrt(N),N,D);    % initial conditions drawn from a normal distribution
%x0 = perfectLactice(N, LinkNumber, D, true, true, (floor(nthroot(N,D)+1))^D); % initial conditions on a correct lattice
%x0 = perfectLactice(N, LinkNumber, D) + randCircle(N, delta, D); % initial conditions on a deformed lattice
%x0 = perfectLactice(N, LinkNumber, D, true, true, (floor(nthroot(N,D)+1))^D ) + randCircle(N, delta, D); % initial conditions on a deformed lattice

speeds0 = abs(normrnd(median(identification.mean_s),median(identification.std_s),N,1));
theta0 = 2*pi*rand(N,1)-pi;
v0 = speeds0 .* [cos(theta0), sin(theta0)];
%v0 = zeros(size(x0));

%% Run Simulation
[xVec, uVec, vVec] = Simulator(x0, v0, Simulation, Dynamics, GlobalIntFunction, LocalIntFunction, Environment);

%% Analysis
if smoothing
    xVec = movmean(xVec,3);
    %xVec = movmean(xVec,3);
end

timeInstants = 0:Simulation.deltaT:Simulation.Tmax;

% derivate quantities
[~, vVec_grad] = gradient(xVec, 1, Simulation.deltaT, 1);
vVec_diff = diff(xVec)/Simulation.deltaT;
speed = vecnorm(vVec,2,3);
speed_grad = vecnorm(vVec_grad,2,3);
speed_diff = [speeds0'; vecnorm(vVec_diff,2,3)];

theta = atan2(vVec(:,:,2), vVec(:,:,1));
for i=1:length(timeInstants)-1
    % angular velocity
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
    links(i)=size(B,2);                                     % number of links of the swarm
    M = buildRigidityMatrix(x, B);                          % rigidity matrix
    
    rigidity(i) = rank(M)==D*N-D*(D+1)/2;                   % check infinitesimal rigidity
end

%% PLOTS

% create output folder, save data and parameters

if outputDir
    counter=1;
    while exist(fullfile(outputDir,[datestr(now, 'yyyy_mm_dd_'),Dynamics.model,'_',num2str(counter)]),'dir')
        counter=counter+1;
    end
    path=fullfile(outputDir, [datestr(now, 'yyyy_mm_dd_'),Dynamics.model,'_',num2str(counter)]);
    mkdir(path)
    disp('Saving data in ' + string(path))
    save(fullfile(path, 'data'))
    
    fileID = fopen(fullfile(path, 'parameters.txt'),'wt');
    fprintf(fileID,'StabilityAnalysis\n\n');
    fprintf(fileID,'Date: %s\n',datestr(now, 'dd/mm/yy'));
    fprintf(fileID,'Time: %s\n\n',datestr(now, 'HH:MM'));
    fprintf(fileID,'Parameters:\n\n');
    fprintf(fileID,'N= %d\n',N);
    fprintf(fileID,'D= %d\n',D);
    fprintStruct(fileID,Simulation)
    fprintf(fileID,'Dynamics:\n');
    fprintStruct(fileID,Dynamics)
    fprintf(fileID,'GlobalIntFunction:\n');
    fprintStruct(fileID,GlobalIntFunction)
    fprintf(fileID,'LocalIntFunction:\n');
    fprintStruct(fileID,LocalIntFunction)
    fprintf(fileID,'smoothing= %s\n',mat2str(smoothing));
    fprintf(fileID,'delta= %.2f\n',delta);
    fclose(fileID);
end

% SWARM
figure
if isfield(LocalIntFunction, 'DistanceRange')
    plotSwarmInit(x0, 0, LocalIntFunction.DistanceRange(1), LocalIntFunction.DistanceRange(2), Simulation.arena);
else
    plotSwarmInit(x0, 0, inf, inf, Simulation.arena);
end
if Simulation.drawTraj; plotTrajectory(xVec, false, [0,0.7,0.9], Simulation.drawTraj); end
if outputDir
    saveas(gcf, fullfile(path, 'trajectories'))
    saveas(gcf, fullfile(path, 'trajectories'),'png')
end

if ~strcmp(GlobalIntFunction.function,'None') % GLOBAL INTERACTION FUNCTION
    figure
    set(gcf,'Position',[100 500 560 420*0.6])
    hold on
    fplot(@(x) globalInteractionForce(x, GlobalIntFunction),[0, 2], 'LineWidth', 1.5)
    plot([1], [0], 'r.','MarkerSize', 25)
    yticks([-1:2])
    xticks(sort([0:0.5:3,Rmax]))
    ylim([-0.6, 2])
    grid on
    set(gca,'FontSize',14)
    ylabel('$f(z)$', 'Interpreter','latex','FontSize',22, 'rotation',0,'VerticalAlignment','middle')
    xlabel('$z$', 'Interpreter','latex','FontSize',22)    
    box
    if outputDir
    saveas(gcf,fullfile(path, 'radial_inter_func'))
    saveas(gcf,fullfile(path, 'radial_inter_func'),'png')
    end
end

if ~strcmp(LocalIntFunction.function, 'None') % LOCAL INTERACTION FUNCTION
    figure 
    hold on
    fplot(@(alfa) localInteractionForce(zeros(1,D), [cos(alfa),sin(alfa)], LocalIntFunction),[-pi/LinkNumber, pi/LinkNumber])
    plot([0], [0], 'r.','MarkerSize', 30)
    ylim([-1.2 1.2])
    xlim([-pi/LinkNumber pi/LinkNumber])
    yticks([-1 0 1])
    xticks([-pi/LinkNumber, 0, pi/LinkNumber])
    set(gca,'XTickLabel',{'-\pi/4','0','\pi/4'})
    grid on
    title('f_n(\theta)')
    xlabel('\theta')
    set(gca,'FontSize',14)
end

figure % TIME PLOT - SPEED and ANGULAR VELOCITY
subplot(2,4,[1 2 3])
plotWithShade(timeInstants, median(speed,2), min(speed, [], 2), max(speed, [], 2), 'b', 0.3);
if isfield(Environment,'EnvUniform')
    highlightInputs(Environment.EnvUniform.Times, Environment.EnvUniform.Values, 'r', 0.25)
end
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
%plotWithShade(timeInstants, median(omega,2), min(omega, [], 2), max(omega, [], 2), 'b', 0.3);
if isfield(Environment,'EnvUniform')
    highlightInputs(Environment.EnvUniform.Times, Environment.EnvUniform.Values, 'r', 0.25)
end
xlabel('t [s]')
ylabel('ang. vel. [rad/s]')
rng=ylim;
box on
subplot(2,4,8)
h=histogram(abs(omega(:)),'Orientation','horizontal');
%h=histogram(omega(:),'Orientation','horizontal');
ylim(rng);
set(gca,'xtick',[])
if outputDir
    saveas(gcf,fullfile(path, 'time_plot'))
    saveas(gcf,fullfile(path, 'time_plot'),'png')
end

figure % SCATTER PLOT - SPEED and ANGULAR VELOCITY
s=scatterhist(speed(:),abs(omega(:)), 'Location','NorthEast','Direction','out');
xlabel(s,'speed')
ylabel(s,'ang. vel. [rad/s]')
s(1).YAxisLocation = 'left';
s(1).XAxisLocation = 'bottom';
s(2).Position = [0.1    0.82   0.7    0.125];
s(3).Position = [0.82   0.1    0.125    0.7];
s(1).Position(3) = 0.7;
s(1).Position(4) = 0.7;
if outputDir
saveas(gcf,fullfile(path, 'scatter_plot'))
saveas(gcf,fullfile(path, 'scatter_plot'),'png')
end

figure % SCATTER PLOT - MEAN SPEED and ANGULAR VELOCITY
s=scatterhist(mean(speed,1),mean(abs(omega),1), 'Location','NorthEast','Direction','out');
xlabel(s,'mean speed [px/s]')
ylabel(s,'mean ang. vel. [rad/s]')
s(1).YAxisLocation = 'left';
s(1).XAxisLocation = 'bottom';
s(2).Position = [0.1    0.82   0.7    0.125];
s(3).Position = [0.82   0.1    0.125    0.7];
s(1).Position(3) = 0.7;
s(1).Position(4) = 0.7;
if outputDir
saveas(gcf,fullfile(path, 'scatter_plot_mean'))
saveas(gcf,fullfile(path, 'scatter_plot_mean'),'png')
end

% figure % CORRELETION PLOT - SPEED and ANGULAR VELOCITY
% corrplot([speed(1:end-1,1),speed(2:end,1),omega(1:end-1,1),omega(2:end,1)],VarNames={"v_k", "v_{k+1}", "\omega_k", "\omega_{k+1}"})
% if outputDir
% saveas(gcf,fullfile(path, 'scatter_plot'))
% saveas(gcf,fullfile(path, 'scatter_plot'),'png')
% end
 
% figure % e_d_max
% set(gca,'FontSize',14)
% set(0, 'DefaultFigureRenderer', 'painters');
% set(gcf,'Position',[100 100 560 420*0.6])
% hold on
% line=plot(timeInstants, e_d_max, 'b');
% yline(Rmax-1,'--','LineWidth',2)
% yticks(sort([0:0.1:1, Rmax-1]))
% ylabel('$e$', 'Interpreter','latex','FontSize',22, 'rotation',0,'VerticalAlignment','middle')
% xlabel('t', 'Interpreter','latex','FontSize',22)
% box
% grid
% if outputDir
%     saveas(gcf,fullfile(path, 'e_d_max'))
%     saveas(gcf,fullfile(path, 'e_d_max'),'png')
% end
% 
% figure % links
% plot(timeInstants,links)
% title('links', 'Interpreter','latex','FontSize',22)
% xlabel('t', 'Interpreter','latex','FontSize',22)
% set(gca,'FontSize',14)
% box
% grid
% if outputDir
%     saveas(gcf,fullfile(path, 'links'))
%     saveas(gcf,fullfile(path, 'links'),'png')
% end
% 
% figure % rigidity
% set(gca,'FontSize',14)
% set(gcf,'Position',[100 100 560 420*0.6])
% hold on
% plot(timeInstants,rigidity,'r')
% axis([-inf inf -0.05 1.05])
% title('$\rho$', 'Interpreter','latex','FontSize',22)
% xlabel('t', 'Interpreter','latex','FontSize',22)
% box
% grid
% if outputDir
%     saveas(gcf,fullfile(path, 'rigidity'))
%     saveas(gcf,fullfile(path, 'rigidity'),'png')
% end


