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

N=1500;

%data_folder = '/Volumes/DOMEPEN/Experiments/2023_06_15_Euglena_1/tracking_2023_10_12';  % off
%data_folder = '/Volumes/DOMEPEN/Experiments/2023_06_15_Euglena_7/tracking_2023_10_16';  % switch10s
%data_folder = '/Volumes/DOMEPEN/Experiments/2023_06_26_Euglena_37/tracking_2023_10_12'; % circle light
%data_folder = '/Volumes/DOMEPEN/Experiments/2023_07_10_Euglena_26/tracking_2024_01_30'; % circle light high denisty
%data_folder = '/Volumes/DOMEPEN/Experiments/2023_07_10_Euglena_21/tracking_2024_01_30'; % circle dark
%data_folder = '/Volumes/DOMEPEN/Experiments/2023_06_26_Euglena_33/tracking_2023_10_12'; % gradient central light
data_folder = '/Volumes/DOMEPEN/Experiments/2023_06_26_Euglena_34/tracking_2023_10_12'; % gradient central dark

id_folder = '/Volumes/DOMEPEN/Experiments/2023_06_15_Euglena_7/tracking_2023_10_16'; % folder with identification data
identification_file_name = 'identification_OLS_ds3_sign.txt';

outputDir = '/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations';

%% Loads experiment data
Simulation.arena = [1920,1080]; % size of the simulation window
Simulation.deltaT = 0.5;        % sampling time step
Simulation.dT =     0.01;   % integration time step
Simulation.Tmax = 180;          % maximum simulation time
timeInstants = [0:Simulation.deltaT:Simulation.Tmax];

% load identification data and instantiate simulated agents
identification=readtable(fullfile(id_folder,identification_file_name));
ids=randsample(length(identification.agents),N, true, ones(length(identification.agents),1));
agents = identification(ids,:); 
Dynamics=struct('model','IndependentSDEsWithInput', ...
    'avgSpeed',agents.mu_s, 'rateSpeed', agents.theta_s, 'sigmaSpeed', agents.sigma_s, 'gainSpeed', agents.alpha_s, 'gainDerSpeed', agents.beta_s,...
    'rateOmega', agents.theta_w, 'sigmaOmega', agents.sigma_w, 'gainOmega', agents.alpha_w, 'gainDerOmega', agents.beta_w,...
    'omega', normrnd(0,agents.std_w,N,1), 'oldInput', zeros(N,1));

% Dynamics=struct('model','IndependentSDEsWithInput', ...
%     'avgSpeed',agents.mu_s, 'rateSpeed', agents.theta_s, 'sigmaSpeed', agents.sigma_s, 'gainSpeed', 0, 'gainDerSpeed', -50,...
%     'rateOmega', agents.theta_w, 'sigmaOmega', agents.sigma_w, 'gainOmega', 0, 'gainDerOmega', 5,...
%     'omega', normrnd(0,agents.std_w,N,1), 'oldInput', zeros(N,1));

% load inputs data
if isfile(fullfile(data_folder,'inputs.txt'))   % time varying inputs
    inputs=load(fullfile(data_folder,'inputs.txt'));
    u=inputs(:,1)/255;              %select blue channel and scale in [0,1]
    Environment.Inputs.Times  = timeInstants;
    Environment.Inputs.Values = u;
else                                            % spatial inputs
    inputs=imread(fullfile(fileparts(data_folder),'patterns_cam/pattern_10.0.jpeg'));
    u=double(inputs(:,:,3))'/255;    %select blue channel and scale in [0,1]
    K=ones(51)/51^2;
    u=conv2(u,K,'same');
    Environment.Inputs.Points = {linspace(-Simulation.arena(1),Simulation.arena(1),size(inputs,2))/2, linspace(-Simulation.arena(2),Simulation.arena(2),size(inputs,1))/2};
    %Environment.Inputs.Values = linspace(-1,1,length(Environment.Inputs.Points{1}))' * ones(1,length(Environment.Inputs.Points{2}));
    Environment.Inputs.Values = u;
end

%% Create Initial Conditions
%rng(1,'twister'); % set the randomn seed to have reproducible results

x0=randCircle(N, 1000, D);                 % initial conditions drawn from a uniform disc

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

xFinal_inWindow = squeeze(xVec(end,(xVec(end,:,1)>-Simulation.arena(1)/2 & xVec(end,:,1)<Simulation.arena(1)/2 ...
                        & xVec(end,:,2)>-Simulation.arena(2)/2 & xVec(end,:,2)<Simulation.arena(2)/2),:));

                    
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
    fprintf(fileID,'SimulateDOMEexp\n\n');
    fprintf(fileID,'Experiment: %s\n',data_folder);
    fprintf(fileID,'Parameters: %s\n\n',identification_file_name);
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
    plotSwarmInit(x0, 0, LocalIntFunction.DistanceRange(1), LocalIntFunction.DistanceRange(2), Simulation.arena);
else
    plotSwarmInit(x0, 0, inf, inf, Simulation.arena);
end
if Simulation.drawTraj; plotTrajectory(xVec, false, [0,0.7,0.9], Simulation.drawTraj); end
if outputDir
    saveas(gcf, fullfile(path, 'x_0'))
    saveas(gcf, fullfile(path, 'x_0'),'png')
end

% SWARM final
figure
if isfield(Environment,'Inputs') && isfield(Environment.Inputs,'Points')
    plotEnvField(Environment.Inputs.Points, Environment.Inputs.Values, Simulation.arena)
end
if isfield(LocalIntFunction, 'DistanceRange')
    plotSwarmInit(xFinal_inWindow, Simulation.Tmax, LocalIntFunction.DistanceRange(1), LocalIntFunction.DistanceRange(2), Simulation.arena);
else
    plotSwarmInit(xFinal_inWindow, Simulation.Tmax, inf, inf, Simulation.arena);
end
if Simulation.drawTraj; plotTrajectory(xVec, false, [0,0.7,0.9], Simulation.drawTraj); end
if outputDir
    saveas(gcf, fullfile(path, 'x_final'))
    saveas(gcf, fullfile(path, 'x_final'),'png')
end

% figure % colored trajectories
% hold on
% colors = get(gca, 'ColorOrder');
% final=Simulation.Tmax;
% window = [-1920/3*2,1920/3,-1080/2,1080/2];
% for i=1:N
%     if xVec(final,i,1) > window(1) && xVec(final,i,1) < window(2) && xVec(final,i,2) > window(3) && xVec(final,i,2) < window(4)
%         c = colors(mod(i-1,7)+1,:);
%         plot(xVec(1:final,i,1),xVec(1:final,i,2), 'color', c); 
%         plot(xVec(final,i,1),xVec(final,i,2),'o', 'color', c, 'MarkerFaceColor', c); 
%     end
% end
% xticks([])
% yticks([])
% axis('equal')
% axis(window)
% box on

% light distribution
if isfield(Environment,'Inputs') && isfield(Environment.Inputs,'Points')
    window = [-Simulation.arena(1),Simulation.arena(1),-Simulation.arena(2),Simulation.arena(2)]/2;
    [density_by_input, bins, norm_slope, c_coeff, coefficents] = agentsDensityByInput(Environment.Inputs.Points, Environment.Inputs.Values, xFinal_inWindow, window);
    
    figure % light distribution
    bar((bins(1:end-1)+bins(2:end))/2,density_by_input, 1)
    hold on
    plot(bins,coefficents(1)+coefficents(2)*bins);
    xlabel('Input intensity')
    ylabel('Density of agents')
    title('Final distribution w.r.t. light intensity')
    yticklabels([]);
    %text(max(bins),max(density_by_input),['\rho=',num2str(c_coeff,3)],'HorizontalAlignment','right','FontSize',14)
    text(max(bins),max(density_by_input),['norm slope=',num2str(norm_slope,3)],'HorizontalAlignment','right','FontSize',14)
    if outputDir
    saveas(gcf, fullfile(path, 'light_distribution'))
    saveas(gcf, fullfile(path, 'light_distribution'),'png')
    end
    
     
    img=imread(fullfile(fileparts(data_folder),'images/fig_180.0.jpeg'));
    img_grey = double(img(:,:,1))/255;
    th=0.3;
    mask=img_grey>th;
    x_vec = linspace(window(1),window(2),size(mask,2));
    y_vec = linspace(window(3),window(4),size(mask,1));
    figure; imagesc(x_vec,y_vec,img_grey); axis('equal')
    figure; imagesc(x_vec,y_vec,mask); axis('equal'); title(['mask=',num2str(th)])
    
    figure % light distribution
    [density_by_input, bins, norm_slope, c_coeff, coefficents] = agentsDensityByInput(Environment.Inputs.Points, Environment.Inputs.Values, mask, window);
    bar((bins(1:end-1)+bins(2:end))/2,density_by_input, 1)
    hold on
    plot(bins,coefficents(1)+coefficents(2)*bins);
    xlabel('Input intensity')
    ylabel('Density of agents')
    title('Final distribution w.r.t. light intensity')
    yticklabels([]);
    %text(max(bins),max(density_by_input),['\rho=',num2str(c_coeff,3)],'HorizontalAlignment','right','FontSize',14)
    text(max(bins),max(density_by_input),['norm slope=',num2str(norm_slope,3)],'HorizontalAlignment','right','FontSize',14)
end

% COMPARE RESULTS
% [MSE_speed,MSE_omega,NMSE_speed,NMSE_omega,NMSE_total] = compareResults({data_folder,path}, path);

% figure % TIME PLOT - SPEED and ANGULAR VELOCITY
% subplot(2,4,[1 2 3])
% plotWithShade(timeInstants, median(speed,2), min(speed, [], 2), max(speed, [], 2), 'b', 0.3);
% if isfield(Environment,'Inputs')
%     highlightInputs(Environment.Inputs.Times, Environment.Inputs.Values, 'r', 0.25)
% end
% xlabel('t [s]')
% ylabel('speed')
% rng=ylim;
% box on
% subplot(2,4,4)
% h=histogram(speed(:),'Orientation','horizontal');
% ylim(rng);
% set(gca,'xtick',[])
% subplot(2,4,[5 6 7])
% plotWithShade(timeInstants, median(abs(omega),2), min(abs(omega), [], 2), max(abs(omega), [], 2), 'b', 0.3);
% %plotWithShade(timeInstants, median(omega,2), min(omega, [], 2), max(omega, [], 2), 'b', 0.3);
% if isfield(Environment,'Inputs')
%     highlightInputs(Environment.Inputs.Times, Environment.Inputs.Values, 'r', 0.25)
% end
% xlabel('t [s]')
% ylabel('ang. vel. [rad/s]')
% rng=ylim;
% box on
% subplot(2,4,8)
% h=histogram(abs(omega(:)),'Orientation','horizontal');
% %h=histogram(omega(:),'Orientation','horizontal');
% ylim(rng);
% set(gca,'xtick',[])
% if outputDir
%     saveas(gcf,fullfile(path, 'time_plot'))
%     saveas(gcf,fullfile(path, 'time_plot'),'png')
% end
% 
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
% saveas(gcf,fullfile(path, 'scatter_plot_mean'))
% saveas(gcf,fullfile(path, 'scatter_plot_mean'),'png')
% end

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


