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
%close all
%clc

%% Parameters

Ntimes=8;              % How many simulations are launched for each configuration

D=2;                    % number of dimensions [2 or 3]

defaultParam;           % load default parameters

seed=-1;                 % seed for random generator, if negative it is not set

%% Loads DOME experiment data
%data_folder = '/Volumes/DOMEPEN/Experiments/2023_06_26_Euglena_33/tracking_2023_10_12'; % gradient central light
data_folder = '/Volumes/DOMEPEN/Experiments/2023_06_26_Euglena_34/tracking_2023_10_12'; % gradient central dark
id_folder = '/Volumes/DOMEPEN/Experiments/2023_06_15_Euglena_7/tracking_2023_10_16'; % folder with identification data
identification_file_name = 'identification_OLS_ds2_sign.txt';
outputDir = '/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/DOME/simulations';

Simulation.arena = [1920,1080]; % size of the simulation window
Simulation.deltaT = 0.5;        % sampling time step
Simulation.dT =     0.01;   % integration time step
Simulation.Tmax = 180;          % maximum simulation time
timeInstants = [0:Simulation.deltaT:Simulation.Tmax];

% load identification data and instantiate simulated agents
identification=readtable(fullfile(id_folder,identification_file_name));
ids=randsample(length(identification.agents),N, true, ones(length(identification.agents),1));
agents = identification(ids,:); 
% Dynamics=struct('model','IndependentSDEsWithInput', ...
%     'avgSpeed',agents.mu_s, 'rateSpeed', agents.theta_s, 'sigmaSpeed', agents.sigma_s, 'gainSpeed', agents.alpha_s, 'gainDerSpeed', agents.beta_s,...
%     'rateOmega', agents.theta_w, 'sigmaOmega', agents.sigma_w, 'gainOmega', agents.alpha_w, 'gainDerOmega', agents.beta_w,...
%     'omega', normrnd(0,agents.std_w,N,1), 'oldInput', zeros(N,1));

Dynamics=struct('model','IndependentSDEsWithInput', ...
    'avgSpeed',agents.mu_s, 'rateSpeed', agents.theta_s, 'sigmaSpeed', agents.sigma_s, 'gainSpeed', 0, 'gainDerSpeed', 0,...
    'rateOmega', agents.theta_w, 'sigmaOmega', agents.sigma_w, 'gainOmega', 0, 'gainDerOmega', 0,...
    'omega', normrnd(0,agents.std_w,N,1), 'oldInput', zeros(N,1));

% load inputs data
if isfile(fullfile(data_folder,'inputs.txt'))   % time varying inputs
    inputs=load(fullfile(data_folder,'inputs.txt'));
    u=inputs(:,1)/255;              %select blue channel and scale in [0,1]
    Environment.Inputs.Times  = timeInstants;
    Environment.Inputs.Values = u;
else                                            % spatial inputs
    inputs=imread(fullfile(fileparts(data_folder),'patterns_cam/pattern_10.0.jpeg'));
    u=double(inputs(:,:,3))'/255;    %select blue channel and scale in [0,1]
    Environment.Inputs.Points = {linspace(-Simulation.arena(1),Simulation.arena(1),size(inputs,2))/2, linspace(-Simulation.arena(2),Simulation.arena(2),size(inputs,1))/2};
    %Environment.Inputs.Values = linspace(-1,1,length(Environment.Inputs.Points{1}))' * ones(1,length(Environment.Inputs.Points{2}));
    Environment.Inputs.Values = u;
end

%% variable parameters
% One or multiple parameters can be modified at the same time.
% Parameters must be existing variables.
% The values specified here overwrite the default ones.
% parameters(1).name='delta';
% parameters(1).values=[0:0.5:1];
% parameters(2).name='N';
% parameters(2).values=[50, 100, 150];

parameters(1).name='Dynamics.gainDerSpeed';
parameters(1).values=[-10,-5,-2,-1,0,1,2,5,10]*5;
%parameters(1).values=[-10,0,10];
parameters(2).name='Dynamics.gainDerOmega';
parameters(2).values=[-1,-0.5,-0.2,-0.1,0,0.1,0.2,0.5,1]*5;
%parameters(2).values=[-10,10];

%% Preallocate
p=cartesianProduct({parameters.values});

Nparameters=length(parameters);
Nconfig=size(p, 1);

timeInstants = 0:Simulation.deltaT:Simulation.Tmax;

xVec=nan(length(timeInstants),N,D);
x_f = nan(Nconfig,Ntimes,N,D);
e_L=nan(Nconfig,Ntimes);
e_theta=nan(Nconfig,Ntimes);
Tr_vec=nan(Nconfig,Ntimes);
success_vec = nan(Nconfig,Ntimes);
e_d_max_vec = nan(Nconfig, Ntimes);
V_vec = nan(Nconfig, Ntimes);
rigid_vec = nan(Nconfig, Ntimes);
norm_slope = nan(Nconfig, Ntimes);



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
        
        evalin('base', [parameters(j).name, '=', num2str(p(i_times,j)), ';'] );
        disp(['> ',parameters(j).name,' = ', num2str(p(i_times,j)) ])
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
        x0Data(k_times,:,:) = randCircle(N, 1000, D);                          % initial conditions drawn from a uniform disc
        %x0Data(k_times,:,:) = normrnd(0,0.1*sqrt(N),N,D);                   % initial conditions drawn from a normal distribution
        %x0Data(k_times,:,:) = perfectLactice(N, LinkNumber, D, true, true, (floor(nthroot(N,D)+1))^D );        % initial conditions on a correct lattice
        %x0Data(k_times,:,:) = perfectLactice(N, LinkNumber, D, true, true, (floor(nthroot(N,D)+1))^D ) + randCircle(N, delta, D); % initial conditions on a deformed lattice
    end
    
    parfor k_times=1:Ntimes
        % run simulation
        [xVec] = Simulator(squeeze(x0Data(k_times,:,:)), v0, Simulation, Dynamics, GlobalIntFunction, LocalIntFunction, Environment);
        
        % analyse final configuration
        xFinal=squeeze(xVec(end,:,:));
        x_f(i_times,k_times,:,:) = xFinal;
        xFinal_inWindow = squeeze(xVec(end,(xVec(end,:,1)>-Simulation.arena(1)/2 & xVec(end,:,1)<Simulation.arena(1)/2 ...
                        & xVec(end,:,2)>-Simulation.arena(2)/2 & xVec(end,:,2)<Simulation.arena(2)/2),:));
        
        B = buildIncidenceMatrix(xFinal, Rmax);
        links(i_times,k_times)=size(B,2);
        M = buildRigidityMatrix(xFinal, B);
        rigid_vec(i_times,k_times) = rank(M)==D*N-D*(D+1)/2;
        e_d_max_vec(i_times,k_times) = getMaxLinkLengthError(xFinal, 1, 0, Rmax);   % max distance from the deisred link length.
        
        F = griddedInterpolant(Environment.Inputs.Points,Environment.Inputs.Values, 'linear', 'nearest');
        envInput = F(xFinal_inWindow(:,1),xFinal_inWindow(:,2));    % input intensity measured by the agents
        x_vec = linspace(-Simulation.arena(1)/2,Simulation.arena(1)/2,100);
        y_vec = linspace(-Simulation.arena(2)/2,Simulation.arena(2)/2,100);
        [x_mesh, y_mesh] = meshgrid(x_vec, y_vec);
        [pixels_by_input,bins] = histcounts(F(x_mesh',y_mesh'), 5);
        [agents_by_input,bins] = histcounts(envInput, bins);
        density_by_input = agents_by_input./pixels_by_input;        
        c_matrix = corrcoef(density_by_input, bins(1:end-1));
        c_coeff = c_matrix(1,2);
        coefficents = [ones(length(bins(1:end-1)),1),(bins(1:end-1)+bins(2:end))'/2]\density_by_input';
        norm_slope(i_times,k_times) = coefficents(2)/mean(density_by_input);
    end
    fprintf('Elapsed time is %.2f s.\n\n',toc)
end


%% Output in command window

fprintf('\n --- \nNtimes=%d seed=%d \n\nConfig \t|',Ntimes,seed)
for j_times=1:Nparameters
    fprintf('%s\t',string(parameters(j_times).name));
end
fprintf('\n')
for i_times=1:Nconfig
    fprintf(' %d \t|', i_times)
    for j_times=1:Nparameters
        fprintf('%.2f\t\t',string(p(i_times,j_times)));
    end
    fprintf('\n')
end

%% Plots

% create folder, save data and parameters
if outputDir
    counter=1;
    while exist(fullfile(outputDir,[datestr(now, 'yyyy_mm_dd_'),Dynamics.model,'_tuning',num2str(counter)]),'dir')
        counter=counter+1;
    end
    path=fullfile(outputDir, [datestr(now, 'yyyy_mm_dd_'),Dynamics.model,'_tuning',num2str(counter)]);
    mkdir(path)
    disp('Saving data in ' + string(path))
    save(fullfile(path, 'data'))
    
    fileID = fopen(fullfile(path, 'parameters.txt'),'wt');
    fprintf(fileID,'SequentialLauncher\n\n');
    fprintf(fileID,'Date: %s\n',datestr(now, 'dd/mm/yy'));
    fprintf(fileID,'Time: %s\n\n',datestr(now, 'HH:MM'));
    fprintf(fileID,'Ntimes= %d\n\n',Ntimes);
    fprintf(fileID,'Parameters:\n\n');
    fprintf(fileID,'N= %d\n',N);
    fprintf(fileID,'D= %d\n\n',D);
    fprintf(fileID,'Simulation parameters:\n');
    fprintStruct(fileID,Simulation)
    fprintf(fileID,'Changing parameters:\n');
    fprintStruct(fileID,parameters)
    fprintf(fileID,'Dynamics:\n');
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
        saveas(gcf,fullfile(path, 'e_rho'))
        saveas(gcf,fullfile(path, 'e_rho'),'png')
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
        saveas(gcf,fullfile(path, 'norm_slope'))
        saveas(gcf,fullfile(path, 'norm_slope'),'png')
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
            plotSwarmInit(squeeze(x_f(sub2ind([length(parameters(1).values), length(parameters(2).values)], i_x, i_y),1,:,:)), Simulation.Tmax, inf, inf, Simulation.arena);
            xticks([]); yticks([])
            title([parameters(1).name,'=' num2str(parameters(1).values(i_x)),' ', parameters(2).name,'=' num2str(parameters(2).values(i_y))])
            title(['\beta_v=' num2str(parameters(1).values(i_x)),' ','\beta_\omega=' num2str(parameters(2).values(i_y))])
        end
    end
    set(gcf,'Position',[100 500 200*swarms_to_show 300*2])
    if outputDir
        saveas(gcf,fullfile(path, 'x'))
        saveas(gcf,fullfile(path, 'x'),'png')
    end
end





