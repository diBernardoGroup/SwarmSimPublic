%
%LauncherMicroorg Set the parameters and launch a simulation of a swarm of microorganisms.
%
%   See also: Launcher
%
%   Authors:    Andrea Giusti
%   Date:       2023
%

%% Clear environment
close all
clear

%% Parameters

%rng(1,'twister');      % set the randomn seed to have reproducible results
defaultParamMicroorg;   % load default parameters to simulate microorganisms

Render.drawON=true;     % draw swarm during simulation (if N is large slows down the simulation)

%% Run Simulation
[xVec, uVec, vVec] = Simulator(x0, v0, Simulation, Dynamics, Render, GlobalIntFunction, LocalIntFunction, Environment);

%% Analysis
if smoothing
    xVec = movmean(xVec,3);
end

% derivate quantities
% vVec = [xVec(2,:,:)-xVec(1,:,:); diff(xVec)]/Simulation.deltaT;
speed = vecnorm(vVec,2,3);

theta = atan2(vVec(:,:,2), vVec(:,:,1));
omega(1,:) = angleBetweenVectors(squeeze(xVec(2,:,:)-xVec(1,:,:)),squeeze(xVec(3,:,:)-xVec(2,:,:)))';
omega(2,:) = angleBetweenVectors(squeeze(xVec(2,:,:)-xVec(1,:,:)),squeeze(xVec(3,:,:)-xVec(2,:,:)))';
for i=3:length(Simulation.timeInstants)
    % angular velocity
    omega(i,:) = angleBetweenVectors(squeeze(vVec(i-1,:,:)),squeeze(vVec(i,:,:)))';
end
omega=omega/Simulation.deltaT;

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

figure % colored trajectories
hold on
colors = get(gca, 'ColorOrder');
final=ceil(size(xVec,1)/1);
window = [-Simulation.arena(1), Simulation.arena(1), -Simulation.arena(2), Simulation.arena(2)]/2;
for i=1:10:N
    if xVec(final,i,1) > window(1) && xVec(final,i,1) < window(2) && xVec(final,i,2) > window(3) && xVec(final,i,2) < window(4)
        c = colors(mod(i-1,7)+1,:);
        plot(xVec(1:final,i,1),xVec(1:final,i,2), 'color', c); 
        plot(xVec(final,i,1),xVec(final,i,2),'o', 'color', c, 'MarkerFaceColor', c); 
    end
end
xticks([])
yticks([])
axis('equal')
axis(window)
box on
title('Trajectories')

figure % TIME PLOT - SPEED and ANGULAR VELOCITY
subplot(2,4,[1 2 3])
plotWithShade(Simulation.timeInstants, median(speed,2), min(speed, [], 2), max(speed, [], 2), 'b', 0.3);
if isfield(Environment,'Inputs') && isfield(Environment.Inputs,'Times')
    highlightInputs(Environment.Inputs.Times, Environment.Inputs.Values, 'r', 0.25)
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
plotWithShade(Simulation.timeInstants, median(abs(omega),2), min(abs(omega), [], 2), max(abs(omega), [], 2), 'b', 0.3);
%plotWithShade(Simulation.timeInstants, median(omega,2), min(omega, [], 2), max(omega, [], 2), 'b', 0.3);
if isfield(Environment,'Inputs') && isfield(Environment.Inputs,'Times')
    highlightInputs(Environment.Inputs.Times, Environment.Inputs.Values, 'r', 0.25)
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
    saveas(gcf,fullfile(output_path, 'time_plot'))
    saveas(gcf,fullfile(output_path, 'time_plot'),'png')
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
% saveas(gcf,fullfile(output_path, 'scatter_plot'))
% saveas(gcf,fullfile(output_path, 'scatter_plot'),'png')
% end
% 
% figure % SCATTER PLOT - MEAN SPEED and ANGULAR VELOCITY
% s=scatterhist(mean(speed,1),mean(abs(omega),1), 'Location','NorthEast','Direction','out');
% xlabel(s,'mean speed [px/s]')
% ylabel(s,'mean ang. vel. [rad/s]')
% s(1).YAxisLocation = 'left';
% s(1).XAxisLocation = 'bottom';
% s(2).Position = [0.1    0.82   0.7    0.125];
% s(3).Position = [0.82   0.1    0.125    0.7];
% s(1).Position(3) = 0.7;
% s(1).Position(4) = 0.7;
% if outputDir
% saveas(gcf,fullfile(output_path, 'scatter_plot_mean'))
% saveas(gcf,fullfile(output_path, 'scatter_plot_mean'),'png')
% end

% figure % CORRELETION PLOT - SPEED and ANGULAR VELOCITY
% corrplot([speed(1:end-1,1),speed(2:end,1),omega(1:end-1,1),omega(2:end,1)],VarNames={"v_k", "v_{k+1}", "\omega_k", "\omega_{k+1}"})
% if outputDir
% saveas(gcf,fullfile(output_path, 'scatter_plot'))
% saveas(gcf,fullfile(output_path, 'scatter_plot'),'png')
% end

% light distribution
if isfield(Environment,'Inputs') && isfield(Environment.Inputs,'Points')
    F = griddedInterpolant(Environment.Inputs.Points,Environment.Inputs.Values, 'linear', 'nearest');
    envInput = F(xFinal_inWindow(:,1),xFinal_inWindow(:,2));    % input intensity measured by the agents
    
    x_vec = linspace(-Simulation.arena(1)/2,Simulation.arena(1)/2,100);
    y_vec = linspace(-Simulation.arena(2)/2,Simulation.arena(2)/2,100);
    [x_mesh, y_mesh] = meshgrid(x_vec, y_vec);
    [pixels_by_input,bins] = histcounts(F(x_mesh',y_mesh'), 10);
    [agents_by_input,bins] = histcounts(envInput, bins);
    
    figure % light distribution
    bar((bins(1:end-1)+bins(2:end))/2,agents_by_input./pixels_by_input, 1)
    xlabel('Input intensity')
    ylabel('Density of agents')
    title('Final distribution')
        
    figure 
    histogram(envInput,bins)
    xlabel('Input intensity')
    ylabel('Number of agents')
    title('Final distribution')
end

if Render.recordVideo
    copyfile('./Output/video.mp4',output_path)
end
