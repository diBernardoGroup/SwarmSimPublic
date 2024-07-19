%
%SequentialLauncher allows to set multiple values of the parameters and
%   launch multiple simulations, from different intial conditions, for each
%   configuration and analyse the steady state results.
%   If multiple parameters are defined all the combinations are tested.
%   If 1 or 2 parameters are varied the results are plotted.
%   It is used to study the effect of parameters on the system.
%
%   Notes:
%       Running this script can take long time (up to hours)
%       For better performances install the parallel computing toolbox
%
%   See also: Launcher, MultiLauncher
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%

%% Clear environment
clear
close all
clc

%% Settings

Ntimes=1;              % How many simulations are launched for each configuration

defaultParamMicroorg;  % load default parameters

seed=0;                % seed for random generator, if negative it is not set

makeSimFolders = true; % save data of individual simulations

N = 12000;
Environment.boundary = Simulation.arena * 2;

%% Loads DOME experiment data
experiments_folder = '/Volumes/DOMEPEN/Experiments';
experiment = '/comparisons/Euglena_switch_10/combo5';

id_folder = '/Volumes/DOMEPEN/Experiments/comparisons/Euglena_switch_10/combo5';  % folder with identification data
identification_file_name = 'identification_GB_absw_noalpha_narrow.txt';

outputDir = '/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations';

%% variable parameters
% One or multiple parameters can be modified at the same time.
% Parameters must be existing variables.
% The values specified here overwrite the default ones.

% parameters(1).name = 'identification_file_name';
% parameters(1).values = ["identification_GB_absw_noalpha_narrow.txt","identification_manual.txt"];

parameters(1).name = 'experiment';
parameters(1).values = ["2023_06_14_E_6","2023_06_12_E_3","2023_06_14_E_10","2023_06_13_E_16","2023_07_10_E_26","2023_06_13_E_15"];
parameters(1).tags = ["half_half", "grad_centr_light","grad_centr_dark","grad_lateral","circle_light","circle_dark"];

%% Preallocate
p=cartesianProduct({parameters.values});

Nparameters=length(parameters);
Nconfig=size(p, 1);

timeInstants = 0:Simulation.deltaT:Simulation.Tmax;
window = [-Simulation.arena(1),Simulation.arena(1),-Simulation.arena(2),Simulation.arena(2)]/2;

xVec=nan(length(timeInstants),N,D);
x_f = cell(Nconfig,Ntimes);
norm_slope = nan(Nconfig, Ntimes);

if ~isfield(parameters,'tags')
    parameters.tags = string(parameters.values);
end

% make output folder
if outputDir
    counter=1;
    while exist(fullfile(outputDir, [datestr(now, 'yyyy_mm_dd'),'_Tuning_',strrep(char(parameters(1).name),'.','_'),'_',num2str(counter)]),'dir')
        counter=counter+1;
    end
    output_path=fullfile(outputDir, [datestr(now, 'yyyy_mm_dd'),'_Tuning_',strrep(char(parameters(1).name),'.','_'),'_',num2str(counter)]);
    mkdir(output_path)
end

if Nconfig*Ntimes>10 && makeSimFolders
    warning('Numerous configurations, consider setting makeSimFolders to false')
end

%% Simulation
% for each configuration...
for i_times=1:Nconfig
    tic
    disp('Simulations batch ' + string(i_times) + ' of ' + Nconfig + ':')
    
    % assign parameters' value
    for j=1:Nparameters
        args = split(parameters(j).name,'.');
        if length(args) == 1
            assert(exist(parameters(j).name,'var'), ['Parameter ',parameters(j).name,' not present in the workspace'] )
        else
            assert(exist(string(args(1)),'var'), ["Structure "+ string(args(1)) + " not present in the workspace"] )
            assert(isfield(eval(string(args(1))), string(args(2))), ["Structure "+ string(args(1)) + " do not have field " + string(args(2))])
        end
        
        if isa(p(i_times,j),'string')
            evalin('base', strjoin([parameters(j).name, '="', p(i_times,j), '";'],'') );
        else
            evalin('base', [parameters(j).name, '=', num2str(p(i_times,j)), ';'] );
        end
        disp(['> ',parameters(j).name,' = ', num2str(p(i_times,j)) ])
    end
    
    % load identification data and instantiate simulated agents
    identification=readtable(fullfile(id_folder,identification_file_name));
    ids=randsample(length(identification.agents),N, true, ones(length(identification.agents),1));
    agents = identification(ids,:);
    Dynamics=struct('model','PTWwithSignedInput', ...
        'avgSpeed',agents.mu_s, 'rateSpeed', agents.theta_s, 'sigmaSpeed', agents.sigma_s, 'gainSpeed', agents.alpha_s, 'gainDerSpeed', agents.beta_s,...
        'avgOmega',agents.mu_w, 'rateOmega', agents.theta_w, 'sigmaOmega', agents.sigma_w, 'gainOmega', agents.alpha_w, 'gainDerOmega', agents.beta_w,...
        'omega', normrnd(0,agents.std_w,N,1), 'oldInput', zeros(N,1));
    
    % load inputs data
    experiment = strrep(experiment,'_E_','_Euglena_');
    data_folder = fullfile(experiments_folder, experiment);
    if isfile(fullfile(data_folder,'inputs.txt'))   % time varying inputs
        inputs    = load(fullfile(data_folder,'inputs.txt'));
        speed_exp = load(fullfile(data_folder,'speeds_smooth.txt'));
        omega_exp = load(fullfile(data_folder,'ang_vel_smooth.txt'));
        u=inputs(:,1)/255;              %select blue channel and scale in [0,1]
        Environment.Inputs.Times  = timeInstants;
        Environment.Inputs.Values = u;
    else                                            % spatial inputs
        u = loadInputPattern(data_folder, pattern_blurring);
        Environment.Inputs.Points = {linspace(-Simulation.arena(1),Simulation.arena(1),size(u,1))/2, linspace(-Simulation.arena(2),Simulation.arena(2),size(u,2))/2};
        Environment.Inputs.Values = flip(u,2);
    end
    
    % create initial conditions
    if seed>=0
        rng(seed,'twister'); % reproducible results
    end
    x0Data=nan(Ntimes,N,D);
    %v0 = zeros(N,D);
    speeds0 = abs(normrnd(median(identification.mean_s),median(identification.std_s),N,1));
    theta0 = 2*pi*rand(N,1)-pi;
    v0 = speeds0 .* [cos(theta0), sin(theta0)];
    for k_times=  1:Ntimes
        %x0Data(k_times,:,:) = randCircle(N, 1000, D);                        % initial conditions drawn from a uniform disc
        x0Data(k_times,:,:) = randRect(N, Simulation.arena*2, D);            % initial conditions drawn from a rectangle
        %x0Data(k_times,:,:) = normrnd(0,0.1*sqrt(N),N,D);                   % initial conditions drawn from a normal distribution
        %x0Data(k_times,:,:) = perfectLactice(N, LinkNumber, D, true, true, (floor(nthroot(N,D)+1))^D );        % initial conditions on a correct lattice
        %x0Data(k_times,:,:) = perfectLactice(N, LinkNumber, D, true, true, (floor(nthroot(N,D)+1))^D ) + randCircle(N, delta, D); % initial conditions on a deformed lattice
    end
    
    for k_times=1:Ntimes
        % run simulation
        [xVec] = Simulator(squeeze(x0Data(k_times,:,:)), v0, Simulation, Dynamics, Render, GlobalIntFunction, LocalIntFunction, Environment);
        
        % analyse final configuration
        xFinal=squeeze(xVec(end,:,:));
        x_f{i_times,k_times} = xFinal;
        xFinal_inWindow = squeeze(xVec(end,(xVec(end,:,1)>window(1) & xVec(end,:,1)<window(2) & xVec(end,:,2)>window(3) & xVec(end,:,2)<window(4)),:));
        
        
        if isfield(Environment,'Inputs') && isfield(Environment.Inputs,'Points')
            [density_by_input_sim, bins, norm_sl, c_coeff] = agentsDensityByInput(Environment.Inputs.Points, Environment.Inputs.Values, xFinal_inWindow, window);
            norm_slope(i_times,k_times) = norm_sl;
            
            % compare with experimental result
            mask = detectObjects(data_folder, background_sub, brightness_thresh);
            [density_by_input_exp, bins, norm_slope, c_coeff, coefficents, agents_by_input, pixels_by_input] = agentsDensityByInput(Environment.Inputs.Points, Environment.Inputs.Values, mask, window);
            tvd(i_times,k_times) = 0.5 * norm(density_by_input_sim-squeeze(density_by_input_exp),1); % Total Variation Distance
            
        else
            [~, vVec] = gradient(xVec, 1, Simulation.deltaT, 1);
            speed = vecnorm(vVec,2,3);
            theta = atan2(vVec(:,:,2), vVec(:,:,1));
            for i=1:length(timeInstants)-1
                % angular velocity
                omega(i,:) = angleBetweenVectors(squeeze(vVec(i,:,:)),squeeze(vVec(i+1,:,:)))';
            end
            omega(length(timeInstants),:) = angleBetweenVectors(squeeze(vVec(length(timeInstants)-1,:,:)),squeeze(vVec(length(timeInstants),:,:)))';
            omega=omega/Simulation.deltaT;
            
            overlap = min(size(omega,1),size(omega_exp,1));
            %             NMSE_speed(i_times,k_times) = goodnessOfFit(median(speed,2,'omitnan'), median(speed_exp,2,'omitnan'), 'NMSE');
            %             NMSE_omega(i_times,k_times) = goodnessOfFit(median(abs(omega(1:end-1,:)),2,'omitnan'), median(abs(omega_exp),2,'omitnan'), 'NMSE');
            %             NMSE_total(i_times,k_times) = mean([NMSE_speed(i_times,k_times), NMSE_omega(i_times,k_times)]);
            
            wmape_speed(i_times,k_times) = mape(median(speed(1:overlap,:),2,'omitnan'), median(speed_exp(1:overlap,:),2,'omitnan'),'wMAPE');
            wmape_omega(i_times,k_times) = mape(median(abs(omega(1:overlap,:)),2,'omitnan'), median(abs(omega_exp(1:overlap,:)),2,'omitnan'),'wMAPE');
            wmape_total(i_times,k_times) = mean([wmape_speed(i_times,k_times), wmape_omega(i_times,k_times)]);
        end
        
        % save single simulation results
        if makeSimFolders
            %sim_ouput_path = fullfile(output_path,[datestr(now, 'yyyy_mm_dd_'), char(parameters(1).name),'_',char(parameters(1).tags(i_times)),'_',num2str(k_times)]);
            sim_ouput_path = fullfile(output_path,[strrep(char(parameters(1).name),'.','_'),'_',char(parameters(1).tags(i_times)),'_',num2str(k_times)]);
            mkdir(sim_ouput_path);
            save(fullfile(sim_ouput_path, 'data'),'xVec','Simulation', 'Dynamics', 'GlobalIntFunction', 'LocalIntFunction', 'Environment');
            
            % SWARM final
            [~,indices_inWindow] = getInWindow(squeeze(xVec(end,:,:)), Simulation.arena);
            xFinal_inWindow = squeeze(xVec(end,indices_inWindow,:));
            xSemiFinal_inWindow = squeeze(xVec(end-1,indices_inWindow,:));
            figure(1)
            cla
            if isfield(Environment,'Inputs') && isfield(Environment.Inputs,'Points')
                plotEnvField(Environment.Inputs.Points, Environment.Inputs.Values, Simulation.arena)
            end
            if Simulation.drawTraj; plotTrajectory(xVec, false, [0,0.7,0.9], Simulation.drawTraj); end
            if isfield(LocalIntFunction, 'DistanceRange')
                plotSwarmInit(xFinal_inWindow, Simulation.Tmax, LocalIntFunction.DistanceRange(1), LocalIntFunction.DistanceRange(2), Simulation.arena, Simulation.arena, false, false, false, Simulation.agentShape, Simulation.agentSize, xSemiFinal_inWindow);
            else
                plotSwarmInit(xFinal_inWindow, Simulation.Tmax, inf, inf, Simulation.arena, Simulation.arena, false, false, false, Simulation.agentShape, Simulation.agentSize, xSemiFinal_inWindow);
            end
            if isfield(Environment,'boundary'); plotBoundary(Environment.boundary); end
            saveas(gcf, fullfile(sim_ouput_path, 'x_final'))
            saveas(gcf, fullfile(sim_ouput_path, 'x_final'),'png')
            
            % difference between light distribution
            figure(2) 
            cla
            hold on
            b_exp = bar((bins(1:end-1)+bins(2:end))/2,density_by_input_exp, 1, FaceColor = 'b', FaceAlpha = 0.5);
            b_sim = bar((bins(1:end-1)+bins(2:end))/2,density_by_input_sim, 1, FaceColor = 'k', FaceAlpha = 0.4);
            legend({'REAL','SIMULATED'},'FontSize',14)
            xlabel('Input intensity','FontSize',14)
            ylabel('Density of agents','FontSize',14)
            yticks([0:0.25:1]);
            text(mean(bins),max(density_by_input_exp)*1.10,['TVD=',num2str(tvd(i_times,k_times),'%.2f')],'HorizontalAlignment','center','FontSize',14)
            ylim([0,max(density_by_input_exp)*1.15])
            xlim([-0.1,1.1])
            xticks(round(bins,2))
            box
                saveas(gcf, fullfile(sim_ouput_path, 'difference_light_distribution'))
                saveas(gcf, fullfile(sim_ouput_path, 'difference_light_distribution'),'png')
            
            
        end
        
    end
    fprintf('Elapsed time is %.2f s.\n\n',toc)
end


%% Output in command window

if isfield(Environment,'Inputs') && isfield(Environment.Inputs,'Points')
    metrics_of_interest = {tvd};
    metrics_color = ['b'];
    metrics_tags = ["TVD"];
else
    metrics_of_interest = {wmape_speed, wmape_total, wmape_omega};
    metrics_color = ['b','k','r'];
    metrics_tags = ["wMAPE_v", "wMAPE_{tot}", "wMAPE_\omega"];
end


fprintf('\n --- \nNtimes=%d seed=%d \n\nConfig \t|',Ntimes,seed)
for j_times=1:Nparameters
    fprintf('%s\t',string(parameters(j_times).name));
end
for j_times=1:length(metrics_tags)
    fprintf('%s\t',string(metrics_tags(j_times)));
end
fprintf('\n')
for i_times=1:Nconfig
    fprintf(' %d \t|', i_times)
    for j_times=1:Nparameters
        if isnumeric(p(i_times,j_times))
            fprintf('%.2f\t\t',p(i_times,j_times));
        else
            fprintf('%s\t\t',string(p(i_times,j_times)));
        end
        for k_times=1:length(metrics_tags)
            fprintf('%.2f\t\t',metrics_of_interest{k_times}(i_times,j_times));
        end
    end
    fprintf('\n')
end

%% Plots

% create folder, save data and parameters
if outputDir
    disp('Saving data in ' + string(output_path))
    save(fullfile(output_path, 'data'))
    
    fileID = fopen(fullfile(output_path, 'parameters.txt'),'wt');
    fprintf(fileID,'SequentialLauncher\n\n');
    fprintf(fileID,'Date: %s\n',datestr(now, 'dd/mm/yy'));
    fprintf(fileID,'Time: %s\n\n',datestr(now, 'HH:MM'));
    fprintf(fileID,'Identification: %s\n\n',fullfile(id_folder,identification_file_name));
    fprintf(fileID,'Ntimes= %d\n\n',Ntimes);
    fprintf(fileID,'Parameters:\n\n');
    fprintf(fileID,'N= %d\n',N);
    fprintf(fileID,'D= %d\n\n',D);
    fprintf(fileID,'Simulation parameters:\n');
    fprintStruct(fileID,Simulation)
    fprintf(fileID,'Changing parameters:\n');
    fprintStruct(fileID,parameters)
    fprintf(fileID,'\nDynamics:\n');
    fprintStruct(fileID,Dynamics)
    fprintf(fileID,'Environment:\n');
    fprintStruct(fileID,Environment)
    fprintf(fileID,'GlobalIntFunction:\n');
    fprintStruct(fileID,GlobalIntFunction)
    fprintf(fileID,'LocalIntFunction:\n');
    fprintStruct(fileID,LocalIntFunction)
    fprintf(fileID,'seed= %d\n',seed);
    fclose(fileID);
end



% plot if Nparameters==1
if Nparameters==1
    if isnumeric(parameters(1).values)
        figure %e_d_max and rigidity
        set(0, 'DefaultFigureRenderer', 'painters');
        subplot(2,1,1)
        hold on
        line=plotWithShade(parameters(1).values, mean(e_d_max_vec,2), min(e_d_max_vec,[],2), max(e_d_max_vec,[],2), 'b', 0.1); %e_d_max_mean(:,1),e_d_max_mean(:,2),e_d_max_mean(:,3), 'b', 0.1);
        yline(Rmax-1,'--','LineWidth',2)
        yticks(sort([0:0.1:1, Rmax-1]))
        xticks(parameters(1).values)
        set(gca,'FontSize',14)
        ylabel('$e$', 'Interpreter','latex','FontSize',22, 'rotation',0,'VerticalAlignment','middle')
        box on
        grid
        
        subplot(2,1,2)
        rigidity_line=plot(parameters(1).values, mean(rigid_vec,2),'Marker','o','Color','r','LineWidth',2,'MarkerSize',6);
        xticks(parameters(1).values)
        yticks([0:0.25:1])
        xlabel(parameters(1).name)
        set(gca,'FontSize',14)
        xlabel('$\delta$', 'Interpreter','latex','FontSize',22)
        ylabel('$\rho$', 'Interpreter','latex','FontSize',22, 'rotation',0,'VerticalAlignment','middle')
        box on
        grid
        if outputDir
            saveas(gcf,fullfile(output_path, 'e_rho'))
            saveas(gcf,fullfile(output_path, 'e_rho'),'png')
        end
        
    else
        figure
        hold on
        for i=1:length(metrics_of_interest)
            plots(i,:)=scatter([1:Nconfig]-(length(metrics_of_interest)-1)*0.1+(i-1)*0.2,metrics_of_interest{i},metrics_color(i));
        end
        xticks([1:Nconfig])
        xticklabels(parameters(1).values)
        xticklabels(parameters(1).tags)
        set(gca, 'TickLabelInterpreter', 'none');
        xlim([0,Nconfig+1])
        ylim([0, max([metrics_of_interest{:}],[],'all')*1.1])
        legend(plots(:,1),metrics_tags)
        set(gca,'FontSize',14)
        box on
        if outputDir
            saveas(gcf,fullfile(output_path, 'metrics'))
            saveas(gcf,fullfile(output_path, 'metrics'),'png')
        end
    end
    
elseif Nparameters==2
    % average over the initial conditions
    links_mean = mean(links,2);
    rigid_mean = mean(rigid_vec,2);
    %     e_d_max_mean = mean(e_d_max_vec,2);
    %     e_d_max_map = reshape(e_d_max_mean, [length(parameters(1).values), length(parameters(2).values)]);
    
    norm_slope_mean = mean(norm_slope,2);
    norm_slope_map = reshape(norm_slope_mean, [length(parameters(1).values), length(parameters(2).values)]);
    
    figure
    [~,lplot]=mysurfc(parameters(1).values, parameters(2).values, norm_slope_map);
    xlabel(parameters(1).name)
    xlabel('$\beta_v$','Interpreter','latex','FontSize',18)
    ylabel(parameters(2).name)
    ylabel('$\alpha_\omega$','Interpreter','latex','FontSize',18)
    ylabel('$\beta_\omega$','Interpreter','latex','FontSize',18)
    title('Photoaccumulation Index')
    hold on
    xlim([-inf, inf])
    ylim([-inf, inf])
    set(gca, 'XTick', parameters(1).values);
    set(gca, 'YTick', parameters(2).values);
    set(gca,'FontSize',14)
    caxis([-1,1])
    if outputDir
        saveas(gcf,fullfile(output_path, 'norm_slope'))
        saveas(gcf,fullfile(output_path, 'norm_slope'),'png')
    end
    
    % SWARM
    figure
    swarms_to_show=min([Nconfig, 6]);
    n_x = length(parameters(1).values);
    n_y = length(parameters(2).values);
    f=tiledlayout(n_y,n_x, 'TileSpacing','tight', 'Padding','tight');
    for i_y=1:n_y
        for i_x=1:n_x
            nexttile(sub2ind([length(parameters(1).values), length(parameters(2).values)], i_x, i_y))
            if isfield(Environment,'Inputs') && isfield(Environment.Inputs,'Points')
                plotEnvField(Environment.Inputs.Points, Environment.Inputs.Values, Simulation.arena)
            end
            plotSwarmInit(squeeze(x_f{sub2ind([length(parameters(1).values), length(parameters(2).values)], i_x, i_y),1}(:,:)), Simulation.Tmax, inf, inf, Simulation.arena);
            xticks([]); yticks([])
            title([parameters(1).name,'=' num2str(parameters(1).values(i_x)),' ', parameters(2).name,'=' num2str(parameters(2).values(i_y))])
            title(['\beta_v=' num2str(parameters(1).values(i_x)),' ','\beta_\omega=' num2str(parameters(2).values(i_y))])
        end
    end
    set(gcf,'Position',[100 500 200*swarms_to_show 300*2])
    if outputDir
        saveas(gcf,fullfile(output_path, 'x_final'))
        saveas(gcf,fullfile(output_path, 'x_final'),'png')
    end
end





