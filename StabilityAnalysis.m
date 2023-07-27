%
%StabilityAnalysis Set the parameters and launch multiple simulations from different initial conditions.
%   Also robustness tests can be run, see AgentsRemoval, NoiseTest and dynamicLattice
%
%   See also: Launcher, SequentialLauncher
%   
%   Authors:    Andrea Giusti
%   Date:       2022
%

%% Clear environment
close all
clear
clc

%% Parameters

Ntimes=10;              % How many simulations are launched for each configuration

defaultParam;           % load default parameters


rng(0,'twister');       % set the randomn seed to have reproducible results


%% Run Simulation
disp(['Running ',num2str(Ntimes),' simulations:'])
for rep=1:Ntimes
    %% Create Initial Conditions    
    %x0=randCircle(N, 2);                               % initial conditions drawn from a uniform disc
    %x0 = normrnd(0,0.1*sqrt(N),N,2);                   % initial conditions drawn from a normal distribution
    x0 = perfectLactice(N, LinkNumber, true, true, (sqrt(N)+1)^2);        % initial conditions on a correct lattice
    
    % initial conditions perturbation
    
    %x0 = perturbInitialConditions(x0, delta, Rmax);    % perturbe initial conditions to obtain a desired value of V
    x0 = x0+randCircle(N, delta);                       % perturbe initial conditions appling a random displacement
    
    %% Run Simulation
    [T_r, success, final_e_theta, final_e_L, final_e_d, finalGRadial, finalGNormal, stopTime, xVec(rep,:,:,:)] = Simulator(x0, LinkNumber, G_radial, G_normal, regularity_thresh, compactness_thresh, Tmax, sigma_actuation, sigma_measure, compassBias, drawON, getMetrics, RadialIntFunction, AgentsRemoval, FaultyAgents, MaxSensingRadius, alpha, beta, dynamicLattice, Rmax);
    
    
    %% ANALYSIS
    time_instants = size(xVec,4);
    for i=1:time_instants % for each time instant...
        x=squeeze(xVec(rep,:,:,i));
        
        e_d(rep,i) = getAvgLinkLengthError(x, 1, 0, Rmax);      % avg distance from the deisred link length
        e_d_max(rep,i) = getMaxLinkLengthError(x, 1, 0, Rmax);   % max distance from the deisred link length.
                                                                 % e_d_max<=(Rmax-1)preserves all the links.
                                                                 % e_d_max(0)<= 2*delta
        
        D = buildIncidenceMatrix(x, Rmax);                  % incidence matrix
        m(rep,i)=size(D,2);                                     % number of links
        M = buildRigidityMatrix(x, D);                      % rigidity matrix
        
        rigidity(rep,i) = rank(M)==2*N-3;                       % check infinitesimal rigidity
    end
end

%% PLOTS
    figure % INITIAL CONDITION
    plotSwarmInit(squeeze(xVec(1,:,:,1)), 0, 0, Rmax)
    
    figure % HALF TIME
    plotSwarmInit(squeeze(xVec(1,:,:,ceil(time_instants))), Tmax/2, 0, Rmax)

    figure % FINAL CONDITION
    plotSwarmInit(squeeze(xVec(1,:,:,time_instants)), Tmax, 0, Rmax)

%     figure % e_d
%     plot(linspace(0,min(stopTime,Tmax,'omitnan'),time_instants),e_d)
%     title('e_d')
%     xlabel('t')
%     set(gca,'FontSize',14)
%     box
%     grid
    
    figure % e_d_max
    set(gca,'FontSize',14)
    set(0, 'DefaultFigureRenderer', 'painters');
    set(gcf,'Position',[100 100 560 420*0.6])
    hold on
    line=plotWithShade(linspace(0,min(stopTime,Tmax,'omitnan'),time_instants), mean(e_d_max), min(e_d_max), max(e_d_max), 'b', 0.2);
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

%     figure % m
%     plot(linspace(0,min(stopTime,Tmax,'omitnan'),time_instants),m)
%     title('m', 'Interpreter','latex','FontSize',22)
%     xlabel('t', 'Interpreter','latex','FontSize',22)
%     set(gca,'FontSize',14)
%     box
%     grid
    
    figure % rigidity
    set(gca,'FontSize',14)
    set(gcf,'Position',[100 100 560 420*0.6])
    hold on
    plot(linspace(0,min(stopTime,Tmax,'omitnan'),time_instants),mean(rigidity))
    axis([-inf inf -0.05 1.05])
    title('$\rho$', 'Interpreter','latex','FontSize',22)
    xlabel('t', 'Interpreter','latex','FontSize',22)
    box
    grid
    
%     figure % RADIAL INTERACTION FUNCTION
%     hold on
%     set(gca,'FontSize',14)
%     fplot(@(x) RadialInteractionForce(x, RadialIntFunction),[0, 2], 'LineWidth', 2)
%     plot([1], [0], 'r.','MarkerSize', 30)
%     yticks([-1:0.5:1])
%     xticks([0:0.5:1, Rmax, 1.5:0.5:3])
%     set(gca,'XTickLabel',{[0:0.5:1], 'R_a', [1.5:0.5:3]})
%     ylim([-0.2 1.2])
%     grid on
%     title('$f(z)$','Interpreter','latex', 'FontSize', 22)
%     xlabel('$z$','Interpreter','latex', 'FontSize', 22)
%     box on
%     grid on

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
