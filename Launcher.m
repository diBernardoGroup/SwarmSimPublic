%
%Launcher Set the parameters and launch a single simulation of the swarm.
%   Also robustness tests can be run, see AgentsRemoval, NoiseTest and dynamicLattice
%
%   See also: LauncherMicroorg, MultiLauncher, SequentialLauncher
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%

%% Clear environment
close all
clear

%% Parameters

%rng(1,'twister');      % set the randomn seed to have reproducible results
defaultParam;           % load default parameters and set initial conditions

Render.drawON=true;     % draw swarm during simulation (if N is large slows down the simulation)

%% Run Simulation
[xVec, uVec, vVec] = Simulator(x0, v0, Simulation, Dynamics, Render, GlobalIntFunction, LocalIntFunction, Environment);

%% Analysis

% metrics
for i=1:length(Simulation.timeInstants) % for each time instant...
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
close all

[~,indices_inWindow] = getInWindow(squeeze(xVec(end,:,:)), Render.window);
xFinal_inWindow = squeeze(xVec(end,indices_inWindow,:));
xSemiFinal_inWindow = squeeze(xVec(end-1,indices_inWindow,:));

% create output folder, save data and parameters
if outputDir
    counter=1;
    while exist(fullfile(outputDir,[datestr(now, 'yyyy_mm_dd_'),Dynamics.model,'_',num2str(counter)]),'dir')
        counter=counter+1;
    end
    output_path=fullfile(outputDir, [datestr(now, 'yyyy_mm_dd_'),Dynamics.model,'_',num2str(counter)]);
    mkdir(output_path)
    disp('Saving data in ' + string(output_path))
    save(fullfile(output_path, 'data'))
    
    fileID = fopen(fullfile(output_path, 'parameters.txt'),'wt');
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

% SWARM initial
figure
if isfield(Environment,'Inputs') && isfield(Environment.Inputs,'Points')
    plotEnvField(Environment.Inputs.Points, Environment.Inputs.Values, Simulation.arena)
end
if isfield(LocalIntFunction, 'DistanceRange')
    plotSwarmInit(x0, 0, LocalIntFunction.DistanceRange(1), LocalIntFunction.DistanceRange(2), Render.window, [Render.window(2)-Render.window(1), Render.window(4)-Render.window(3)]/2, false, false, false, Render.agentShape, Render.agentSize, Render.agentsColor, squeeze(xVec(2,:,:)));
else
    plotSwarmInit(x0, 0, inf, inf, Render.window, [Render.window(2)-Render.window(1), Render.window(4)-Render.window(3)]/2, false, false, false, Render.agentShape, Render.agentSize, Render.agentsColor, squeeze(xVec(2,:,:)));
end
if Render.drawTraj; plotTrajectory(xVec, false, [0,0.7,0.9], Render.drawTraj); end
if outputDir
    saveas(gcf, fullfile(output_path, 'x0'))
    saveas(gcf, fullfile(output_path, 'x0'),'png')
end

% SWARM final
figure
if isfield(Environment,'Inputs') && isfield(Environment.Inputs,'Points')
    plotEnvField(Environment.Inputs.Points, Environment.Inputs.Values, Simulation.arena)
end
if isfield(LocalIntFunction, 'DistanceRange')
    plotSwarmInit(xFinal_inWindow, Simulation.Tmax, LocalIntFunction.DistanceRange(1), LocalIntFunction.DistanceRange(2), Render.window, [Render.window(2)-Render.window(1), Render.window(4)-Render.window(3)]/2, false, false, false, Render.agentShape, Render.agentSize, Render.agentsColor, xSemiFinal_inWindow);
else
    plotSwarmInit(xFinal_inWindow, Simulation.Tmax, inf, inf, Render.window, [Render.window(2)-Render.window(1), Render.window(4)-Render.window(3)]/2, false, false, false, Render.agentShape, Render.agentSize, Render.agentsColor, xSemiFinal_inWindow);
end
if Render.drawTraj; plotTrajectory(xVec, false, [0,0.7,0.9], Render.drawTraj); end
if outputDir
    saveas(gcf, fullfile(output_path, 'xf'))
    saveas(gcf, fullfile(output_path, 'xf'),'png')
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
    title('Interaction function')
    if outputDir
    saveas(gcf,fullfile(output_path, 'radial_inter_func'))
    saveas(gcf,fullfile(output_path, 'radial_inter_func'),'png')
    end
end

% if ~strcmp(LocalIntFunction.function, 'None') % LOCAL INTERACTION FUNCTION
%     figure 
%     hold on
%     fplot(@(alfa) localInteractionForce(zeros(1,D), [cos(alfa),sin(alfa)], LocalIntFunction),[-pi/LinkNumber, pi/LinkNumber])
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
% end
 
figure % e_d_max
set(gca,'FontSize',14)
set(0, 'DefaultFigureRenderer', 'painters');
set(gcf,'Position',[100 100 560 420*0.6])
hold on
line=plot(Simulation.timeInstants, e_d_max, 'b');
yline(Rmax-1,'--','LineWidth',2)
yticks(sort([0:0.1:1, Rmax-1]))
ylabel('$e$', 'Interpreter','latex','FontSize',22, 'rotation',0,'VerticalAlignment','middle')
xlabel('t', 'Interpreter','latex','FontSize',22)
box
grid
if outputDir
    saveas(gcf,fullfile(output_path, 'e_d_max'))
    saveas(gcf,fullfile(output_path, 'e_d_max'),'png')
end

figure % links
plot(Simulation.timeInstants,links)
title('links', 'Interpreter','latex','FontSize',22)
xlabel('t', 'Interpreter','latex','FontSize',22)
set(gca,'FontSize',14)
box
grid
if outputDir
    saveas(gcf,fullfile(output_path, 'links'))
    saveas(gcf,fullfile(output_path, 'links'),'png')
end

figure % rigidity
set(gca,'FontSize',14)
set(gcf,'Position',[100 100 560 420*0.6])
hold on
plot(Simulation.timeInstants,rigidity,'r')
axis([-inf inf -0.05 1.05])
title('$\rho$', 'Interpreter','latex','FontSize',22)
xlabel('t', 'Interpreter','latex','FontSize',22)
box
grid
if outputDir
    saveas(gcf,fullfile(output_path, 'rigidity'))
    saveas(gcf,fullfile(output_path, 'rigidity'),'png')
end

if Render.recordVideo
    copyfile('./Output/video.mp4',output_path)
end
