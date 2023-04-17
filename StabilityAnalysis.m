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

Ntimes=1;              % How many simulations are launched for each configuration

defaultParam;           % load default parameters


rng(0,'twister');       % set the randomn seed to have reproducible results


%% Run Simulation
disp(['Running ',num2str(Ntimes),' simulations:'])
for rep=1:Ntimes
    %% Create Initial Conditions    
    x0=randCircle(N, 2);                               % initial conditions drawn from a uniform disc
    %x0 = normrnd(0,0.1*sqrt(N),N,2);                   % initial conditions drawn from a normal distribution
    %x0 = perfectLactice(N, LinkNumber, true, true, (sqrt(N)+1)^2);        % initial conditions on a correct lattice
    
    %% Run Simulation
    [xVec(rep,:,:,:), vVec(rep,:,:,:), stopTime] = Simulator(x0, Tmax, drawON, getMetrics, GlobalIntFunction, LocalIntFunction);
    
    
    %% ANALYSIS
    time_instants = size(xVec,4);
    for i=1:time_instants % for each time instant...
        x=squeeze(xVec(rep,:,:,i));

    end
end

%% PLOTS
    figure % INITIAL CONDITION
    plotSwarmInit(squeeze(xVec(1,:,:,1)), 0, 0, Rmax)
    
    figure % HALF TIME
    plotSwarmInit(squeeze(xVec(1,:,:,ceil(time_instants))), Tmax/2, 0, Rmax)

    figure % FINAL CONDITION
    plotSwarmInit(squeeze(xVec(1,:,:,time_instants)), Tmax, 0, Rmax)
    
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
