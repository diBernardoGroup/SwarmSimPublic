%
%Launcher Set the parameters and launch a single simulation of the swarm.
%   Also robustness tests can be run, see AgentsRemoval, NoiseTest and dynamicLattice
%
%   See also: BruteForceTuning, SequentialLauncher, StabilityAnalysis
%   
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%

%% Clear environment
close all
clear
clc

%% Parameters

defaultParam;   % load default parameters

LinkNumber=4;   % number of links (6=triangular lattice, 4=square lattice, 3=hexagonal lattice) (L)

RadialIntFunction=struct('function','Lennard-Jones','parameters',[0.15, 5]);

% control gains
G_radial= 15;   % default value for square lattice 15 (G_r)
G_normal = 8;   % default value for square lattice  8 (G_n)


Tmax=30;    % maximum simulation time (simulation is stopped earlier if steady state is reached)

%output options
drawON=true;        % draw swarm during simulation (if N is large slows down the simulation)
getMetrics=true;    % acquire metrics during the simulation (getMetrics=false discard settling times and stop times)

%% Create Initial Conditions
%rng(1,'twister'); % set the randomn seed to have reproducible results

x0=randCircle(N, 2);                 % initial conditions drawn from a uniform disc
%x0 = normrnd(0,0.1*sqrt(N),N,2);    % initial conditions drawn from a normal distribution
%x0 = perfectLactice(N, LinkNumber); % initial conditions on a correct lattice


%% Run Simulation
[T_r, success, final_e_theta, final_e_L, final_e_d, finalGRadial, finalGNormal, stopTime] = Simulator(x0, LinkNumber, G_radial, G_normal, regularity_thresh, compactness_thresh, Tmax, sigma, drawON, getMetrics, RadialIntFunction, AgentsRemoval, NoiseTest, MaxSensingRadius, alpha, beta, dynamicLattice, Rmax);


%% PLOTS
    figure % RADIAL INTERACTION FUNCTION
    hold on
    fplot(@(x) RadialInteractionForce(x, RadialIntFunction),[0, 3])
    plot([1], [0], 'r.','MarkerSize', 40)
    yticks([-1 0 1])
    xticks([0:3])
    ylim([-0.2 1.2])
    grid on
    title('f_r(d)')
    xlabel('d')
    set(gca,'FontSize',14)

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

