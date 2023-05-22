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

%% Parameters

Ntimes=10;              % How many simulations are launched for each configuration

D=2;                    % number of dimensions [2 or 3]

defaultParam;           % load default parameters

seed=0;                 % seed for random generator, if negative it is not set

%% variable parameters
% One or multiple parameters can be modified at the same time.
% Parameters must be existing variables.
% The values specified here overwrite the default ones.
% parameters(1).name='delta';
% parameters(1).values=[0:0.5:1];
% parameters(2).name='N';
% parameters(2).values=[50, 100, 150];

parameters(1).name='GlobalIntFunction.Gain';
parameters(1).values=[0:0.5:2];
parameters(2).name='N';
parameters(2).values=[50, 100, 150];

%% Preallocate
p=cartesianProduct({parameters.values});

Nparameters=length(parameters);
Nconfig=size(p, 1);

timeInstants = 0:Simulation.deltaT:Simulation.Tmax;

xVec=nan(length(timeInstants),N,D);
e_L=nan(Nconfig,Ntimes);
e_theta=nan(Nconfig,Ntimes);
Tr_vec=nan(Nconfig,Ntimes);
success_vec = nan(Nconfig,Ntimes);
e_d_max_vec = nan(Nconfig, Ntimes);
V_vec = nan(Nconfig, Ntimes);
rigid_vec = nan(Nconfig, Ntimes);

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
    v0 = zeros(N,D);
    parfor k_times=  1:Ntimes
        %x0Data(k_times,:,:) = randCircle(N, 2, D);                          % initial conditions drawn from a uniform disc
        %x0Data(k_times,:,:) = normrnd(0,0.1*sqrt(N),N,D);                   % initial conditions drawn from a normal distribution
        %x0Data(k_times,:,:) = perfectLactice(N, LinkNumber, D, true, true, (floor(nthroot(N,D)+1))^D );        % initial conditions on a correct lattice
        x0Data(k_times,:,:) = perfectLactice(N, LinkNumber, D, true, true, (floor(nthroot(N,D)+1))^D ) + randCircle(N, delta, D); % initial conditions on a deformed lattice
    end
    
    parfor k_times=1:Ntimes
        % run simulation
        [xVec] = Simulator(squeeze(x0Data(k_times,:,:)), v0, Simulation, Dynamics, GlobalIntFunction, LocalIntFunction);
        
        % analyse final configuration
        x_f=squeeze(xVec(end,:,:));
        
        B = buildIncidenceMatrix(x_f, Rmax);
        links(i_times,k_times)=size(B,2);
        M = buildRigidityMatrix(x_f, B);
        rigid_vec(i_times,k_times) = rank(M)==D*N-D*(D+1)/2;
        e_d_max_vec(i_times,k_times) = getMaxLinkLengthError(x_f, 1, 0, Rmax);   % max distance from the deisred link length.
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
    while exist(fullfile(outputDir,[datestr(now, 'yyyy_mm_dd_'),Dynamics.model,'_',num2str(counter)]),'dir')
        counter=counter+1;
    end
    path=fullfile(outputDir, [datestr(now, 'yyyy_mm_dd_'),Dynamics.model,'_',num2str(counter)]);
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
    e_d_max_mean = mean(e_d_max_vec,2);
    e_d_max_map = reshape(e_d_max_mean, [length(parameters(1).values), length(parameters(2).values)]);

    figure % e_d_max
    [~,lplot]=mysurfc(parameters(1).values, parameters(2).values, e_d_max_map);
    xlabel(parameters(1).name)
    ylabel(parameters(2).name)
    title('$e$', 'interpreter', 'latex')
    hold on
    xlim([-inf, inf])
    ylim([-inf, inf])
    set(gca, 'XTick', parameters(1).values);
    set(gca, 'YTick', parameters(2).values);
    set(gca,'FontSize',14)
    if outputDir
        saveas(gcf,fullfile(path, 'e'))
        saveas(gcf,fullfile(path, 'e'),'png')
    end
    
end





