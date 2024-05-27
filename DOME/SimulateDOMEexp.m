%
%SimulateDOMEexp Set the parameters and launch a simulation of the swarm.
%
%   See also: Launcher
%
%   Authors:    Andrea Giusti
%   Date:       2024
%

%% Clear environment
close all
clear

%% Parameters

defaultParamMicroorg;               % load default parameters yo simulate microorganisms

% tag='switch_10'; data_folder = '/Volumes/DOMEPEN/Experiments/2023_07_10_Euglena_15/tracking_2023_10_12';  % switch10s
% tag='switch_10'; data_folder = '/Volumes/DOMEPEN/Experiments/comparisons/Euglena_switch_10/combo3';  % switch10s combo
tag='switch_10'; data_folder = '/Volumes/DOMEPEN/Experiments/comparisons/Euglena_switch_10/combo5';  % switch10s combo 5
% tag='switch_5'; data_folder = '/Volumes/DOMEPEN/Experiments/comparisons/Euglena_switch_5/combo';  % switch5s combo
% tag='switch_1'; data_folder = '/Volumes/DOMEPEN/Experiments/comparisons/Euglena_switch_1/combo';  % switch1s combo
% tag='75_ON'; data_folder = '/Volumes/DOMEPEN/Experiments/comparisons/Euglena_75_ON/combo';  % OFF-ON-OFF 75 combo
% tag='150_ON'; data_folder = '/Volumes/DOMEPEN/Experiments/comparisons/Euglena_150_ON/combo';  % OFF-ON-OFF 150 combo
% tag='255_ON'; data_folder = '/Volumes/DOMEPEN/Experiments/comparisons/Euglena_255_ON/combo';  % OFF-ON-OFF 255 combo
% tag='OFF'; data_folder = '/Volumes/DOMEPEN/Experiments/comparisons/Euglena_OFF/combo';  % OFF combo
% tag='ramp'; data_folder = '/Volumes/DOMEPEN/Experiments/comparisons/Euglena_ramp/combo';  % ramp combo

id_folder = '/Volumes/DOMEPEN/Experiments/comparisons/Euglena_switch_10/combo5';  % folder with identification data
%id_folder = '/Volumes/DOMEPEN/Experiments/comparisons/Euglena_switch_10/combo5_old';  % folder with identification data
identification_file_name = 'identification_GA_lim_nomu.txt';

% outputDir = '/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations';
outputDir = '/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations/comparison/Identifications';
%% Loads experiment data

% load identification data and instantiate simulated agents
identification=readtable(fullfile(id_folder,identification_file_name));
ids=randsample(length(identification.agents),N, true, ones(length(identification.agents),1));
agents = identification(ids,:); 
% Dynamics=struct('model','PTWwithInput', ...
%     'avgSpeed',agents.mu_s, 'rateSpeed', agents.theta_s, 'sigmaSpeed', agents.sigma_s, 'gainSpeed', agents.alpha_s, 'gainDerSpeed', agents.beta_s,...
%     'avgOmega',agents.mu_w, 'rateOmega', agents.theta_w, 'sigmaOmega', agents.sigma_w, 'gainOmega', agents.alpha_w, 'gainDerOmega', agents.beta_w,...
%     'omega', normrnd(0,agents.std_w,N,1), 'oldInput', zeros(N,1));

Dynamics=struct('model','PTWwithSignedInput', ...
    'avgSpeed',agents.mu_s, 'rateSpeed', agents.theta_s, 'sigmaSpeed', agents.sigma_s, 'gainSpeed', agents.alpha_s, 'gainDerSpeed', agents.beta_s,...
    'avgOmega',agents.mu_w, 'rateOmega', agents.theta_w, 'sigmaOmega', agents.sigma_w, 'gainOmega', agents.alpha_w, 'gainDerOmega', agents.beta_w,...
    'omega', normrnd(0,agents.std_w,N,1), 'oldInput', zeros(N,1));

% % parameters for strong photodispersion
% Dynamics=struct('model','PTWwithSignedInput', ...
%     'avgSpeed',agents.mu_s, 'rateSpeed', agents.theta_s, 'sigmaSpeed', agents.sigma_s, 'gainSpeed', 0, 'gainDerSpeed', -50,...
%     'rateOmega', agents.theta_w, 'sigmaOmega', agents.sigma_w, 'gainOmega', 0, 'gainDerOmega', 5,...
%     'omega', normrnd(0,agents.std_w,N,1), 'oldInput', zeros(N,1));

% load inputs data
if isfile(fullfile(data_folder,'inputs.txt'))   % time varying inputs
    inputs=load(fullfile(data_folder,'inputs.txt'));
    u=inputs(:,1)/255;              %select blue channel and scale in [0,1]
    Environment.Inputs.Times  = Simulation.timeInstants;
    Environment.Inputs.Values = u;
else                                            % spatial inputs
    [mask, u]= analyseDOMEspatial(fileparts(data_folder), background_sub, brightness_thresh);
    Environment.Inputs.Points = {linspace(-Simulation.arena(1),Simulation.arena(1),size(u,1))/2, linspace(-Simulation.arena(2),Simulation.arena(2),size(u,2))/2};
    Environment.Inputs.Values = flip(u,2);
end

%% Create Initial Conditions
%rng(1,'twister'); % set the randomn seed to have reproducible results

x0=randCircle(N, 1000, D);                 % initial conditions drawn from a uniform disc

speeds0 = abs(normrnd(median(identification.mean_s),median(identification.std_s),N,1));
theta0 = 2*pi*rand(N,1)-pi;
v0 = speeds0 .* [cos(theta0), sin(theta0)];
%v0 = zeros(size(x0));

%% Run Simulation
[xVec, uVec, ~] = Simulator(x0, v0, Simulation, Dynamics, GlobalIntFunction, LocalIntFunction, Environment);

%% Analysis
if smoothing
    xVec = movmean(xVec,3);
    %xVec = movmean(xVec,3);
end

% derivate quantities
[~, vVec_grad] = gradient(xVec, 1, Simulation.deltaT, 1);
vVec_fe = [diff(xVec); xVec(end,:,:)-xVec(end-1,:,:)]/Simulation.deltaT;
vVec_be = [xVec(2,:,:)-xVec(1,:,:); diff(xVec)]/Simulation.deltaT;
speed_grad = vecnorm(vVec_grad,2,3);
speed_fe = vecnorm(vVec_fe,2,3);
speed_be = vecnorm(vVec_be,2,3);
speed = speed_be;

theta = atan2(vVec_grad(:,:,2), vVec_grad(:,:,1));
% angular velocity - gradient
for i=1:length(Simulation.timeInstants)-1
    omega_grad(i,:) = angleBetweenVectors(squeeze(vVec_grad(i,:,:)),squeeze(vVec_grad(i+1,:,:)))';
end
omega_grad(length(Simulation.timeInstants),:) = angleBetweenVectors(squeeze(vVec_grad(length(Simulation.timeInstants)-1,:,:)),squeeze(vVec_grad(length(Simulation.timeInstants),:,:)))';
omega_grad=omega_grad/Simulation.deltaT;

% angular velocity - Forward Euler
for i=1:length(Simulation.timeInstants)-1
    omega_fe(i,:) = angleBetweenVectors(squeeze(vVec_fe(i,:,:)),squeeze(vVec_fe(i+1,:,:)))';
end
omega_fe(length(Simulation.timeInstants),:) = angleBetweenVectors(squeeze(vVec_fe(length(Simulation.timeInstants)-1,:,:)),squeeze(vVec_fe(length(Simulation.timeInstants),:,:)))';
omega_fe=omega_fe/Simulation.deltaT;

% angular velocity - Backward Euler
omega_be(1,:) = angleBetweenVectors(squeeze(xVec(2,:,:)-xVec(1,:,:)),squeeze(xVec(3,:,:)-xVec(2,:,:)))';
omega_be(2,:) = angleBetweenVectors(squeeze(xVec(2,:,:)-xVec(1,:,:)),squeeze(xVec(3,:,:)-xVec(2,:,:)))';
for i=3:length(Simulation.timeInstants)
    omega_be(i,:) = angleBetweenVectors(squeeze(vVec_be(i-1,:,:)),squeeze(vVec_be(i,:,:)))';
end
omega_be=omega_be/Simulation.deltaT;

% angular velocity - Central Derivative
omega_ce(1,:) = angleBetweenVectors(squeeze(xVec(2,:,:)-xVec(1,:,:)),squeeze(xVec(3,:,:)-xVec(2,:,:)))';
for i=2:length(Simulation.timeInstants)-1
    omega_ce(i,:) = angleBetweenVectors(squeeze(xVec(i,:,:)-xVec(i-1,:,:)),squeeze(xVec(i+1,:,:)-xVec(i,:,:)))';
end
omega_ce(length(Simulation.timeInstants),:) = angleBetweenVectors(squeeze(xVec(end-1,:,:)-xVec(end-2,:,:)),squeeze(xVec(end,:,:)-xVec(end-1,:,:)))';
omega_ce=omega_ce/Simulation.deltaT;

omega = omega_be;

xFinal_inWindow = squeeze(xVec(end,(xVec(end,:,1)>-Simulation.arena(1)/2 & xVec(end,:,1)<Simulation.arena(1)/2 ...
                        & xVec(end,:,2)>-Simulation.arena(2)/2 & xVec(end,:,2)<Simulation.arena(2)/2),:));

                    
%% PLOTS

% create output folder, save data and parameters

if outputDir
    counter=1;
    if ~tag 
        tag=Dynamics.model;
    end
    while exist(fullfile(outputDir,[datestr(now, 'yyyy_mm_dd_'),tag,'_',num2str(counter)]),'dir')
        counter=counter+1;
    end
    output_path=fullfile(outputDir, [datestr(now, 'yyyy_mm_dd_'),tag,'_',num2str(counter)]);
    mkdir(output_path)
    disp('Saving data in ' + string(output_path))
    save(fullfile(output_path, 'data'))
    
    fileID = fopen(fullfile(output_path, 'parameters.txt'),'wt');
    fprintf(fileID,'SimulateDOMEexp\n\n');
    fprintf(fileID,'Experiment: %s\n',data_folder);
    fprintf(fileID,'Identification: %s\n\n',fullfile(id_folder,identification_file_name));
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

if isfield(Environment,'Inputs') && isfield(Environment.Inputs,'Points')
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
    saveas(gcf, fullfile(output_path, 'x_0'))
    saveas(gcf, fullfile(output_path, 'x_0'),'png')
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
    saveas(gcf, fullfile(output_path, 'x_final'))
    saveas(gcf, fullfile(output_path, 'x_final'),'png')
end
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

% SPATIAL INPUTS
if isfield(Environment,'Inputs') && isfield(Environment.Inputs,'Points')
    window = [-Simulation.arena(1),Simulation.arena(1),-Simulation.arena(2),Simulation.arena(2)]/2;
    [density_by_input_sim, bins, norm_slope_sim, c_coeff_sim, coefficents, ~,~, u_values_sim] = agentsDensityByInput(Environment.Inputs.Points, Environment.Inputs.Values, xFinal_inWindow, window);
    
    figure % simulation light distribution
    bar((bins(1:end-1)+bins(2:end))/2,density_by_input_sim, 1)
    hold on
    plot(bins,coefficents(1)+coefficents(2)*bins,LineWidth=2);
    xlabel('Input intensity')
    ylabel('Density of agents')
    yticks([0:0.25:1]);
    text(max(bins),max(density_by_input_sim)*1.1,['\rho=',num2str(c_coeff_sim,'%.2f')],'HorizontalAlignment','right','FontSize',14)
    text(max(bins),max(density_by_input_sim)*1.05,['norm slope=',num2str(norm_slope_sim,'%.2f')],'HorizontalAlignment','right','FontSize',14)
    ylim([0,max(density_by_input_sim)*1.15])
    xlim([-0.1,1.1])
    xticks(round(bins,2))
    title('Simulated light distribution')
    if outputDir
        saveas(gcf, fullfile(output_path, 'light_distribution'))
        saveas(gcf, fullfile(output_path, 'light_distribution'),'png')
    end
    
    % get distribution wrt light intensity
    [density_by_input_exp, bins, norm_slope_exp, c_coeff_exp, coefficents, ~,~, u_values_exp] = agentsDensityByInput(Environment.Inputs.Points, Environment.Inputs.Values, mask, window);

    figure
    x_vec = linspace(window(1),window(2),size(mask,2));
    y_vec = linspace(window(3),window(4),size(mask,1));
    box on    
    hold on
    cmap = linspace2([1,1,1], [1,0.5,0.5], 100)';
    colormap(cmap)
    imagesc(x_vec,y_vec,u')
    I=imagesc(x_vec,y_vec,cat(3,zeros(size(mask)),zeros(size(mask)),mask));
    set(I, 'AlphaData', mask);
    axis('equal')
    axis(window)
    xticks([])
    yticks([])
    title('Experimental')
    if outputDir
        saveas(gcf, fullfile(output_path, 'exp_positions'))
        saveas(gcf, fullfile(output_path, 'exp_positions'),'png')
    end
    
    figure % experimental light distribution
    bar((bins(1:end-1)+bins(2:end))/2,density_by_input_exp, 1)
    hold on
    plot(bins,coefficents(1)+coefficents(2)*bins,LineWidth=2);
    xlabel('Input intensity')
    ylabel('Density of agents')
    yticks([0:0.25:1]);
    text(max(bins),max(density_by_input_exp)*1.1,['\rho=',num2str(c_coeff_exp,'%.2f')],'HorizontalAlignment','right','FontSize',14)
    text(max(bins),max(density_by_input_exp)*1.05,['norm slope=',num2str(norm_slope_exp,'%.2f')],'HorizontalAlignment','right','FontSize',14)
    ylim([0,max(density_by_input_exp)*1.15])
    xlim([-0.1,1.1])
    xticks(round(bins,2))
    title('Experimental light distribution')
    box
    if outputDir
        saveas(gcf, fullfile(output_path, 'exp_light_distribution'))
        saveas(gcf, fullfile(output_path, 'exp_light_distribution'),'png')
    end
    
    figure % difference between light distribution
    tvd = 0.5 * norm(density_by_input_exp-density_by_input_sim,1); % Total Variation Distance
    hold on
    b_exp = bar((bins(1:end-1)+bins(2:end))/2,density_by_input_exp, 1, FaceColor = 'b', FaceAlpha = 0.5);
    b_sim = bar((bins(1:end-1)+bins(2:end))/2,density_by_input_sim, 1, FaceColor = 'k', FaceAlpha = 0.4);
    %[f,xi] = ksdensity(u_values_exp, support=[-0.001,1.001], BoundaryCorrection='reflection');
    %f=f/sum(f);
    %plot(xi,f)
    legend({'REAL','SIMULATED'},'FontSize',14)
    xlabel('Input intensity','FontSize',14)
    ylabel('Density of agents','FontSize',14)
    yticks([0:0.25:1]);
    text(mean(bins),max(density_by_input_exp)*1.10,['TVD=',num2str(tvd,'%.2f')],'HorizontalAlignment','center','FontSize',14)
    ylim([0,max(density_by_input_exp)*1.15])
    xlim([-0.1,1.1])
    xticks(round(bins,2))
    box
    if outputDir
        saveas(gcf, fullfile(output_path, 'difference_light_distribution'))
        saveas(gcf, fullfile(output_path, 'difference_light_distribution'),'png')
    end
    
else % TEMPORAL INPUTS

    [metrics_of_interest] = compareResults({data_folder,output_path}, output_path);

end

% figure % TIME PLOT - SPEED and ANGULAR VELOCITY
% subplot(2,4,[1 2 3])
% plotWithShade(Simulation.timeInstants, median(speed,2), min(speed, [], 2), max(speed, [], 2), 'b', 0.3);
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
% plotWithShade(Simulation.timeInstants, median(abs(omega),2), min(abs(omega), [], 2), max(abs(omega), [], 2), 'b', 0.3);
% %plotWithShade(Simulation.timeInstants, median(omega,2), min(omega, [], 2), max(omega, [], 2), 'b', 0.3);
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
%     saveas(gcf,fullfile(output_path, 'time_plot'))
%     saveas(gcf,fullfile(output_path, 'time_plot'),'png')
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

if Simulation.recordVideo
    copyfile('./Output/video.mp4',output_path)
end

