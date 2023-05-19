%
%defaultParam Set the default values of the parameters.
%   It is called by launcher scripts such as Launcher.
%
%   See also: Launcher, SequentialLauncher
%   
%   Authors:    Andrea Giusti
%   Date:       2023
%

%% Default Parameters

outputDir='/Users/andrea/Library/CloudStorage/OneDrive-UniversitàdiNapoliFedericoII/Andrea_Giusti/Projects/stability of geometric lattices/simulations';
%outputDir='';

N=100;                      %number of agents (N)
LinkNumber=6*(D-1);         %number of links (6=triangular lattice, 4=square lattice, 3=hexagonal lattice) (L)


Rmax= (sqrt(5-D)+1)/2;      % maximum lenght of a link (R_a). Must be in [1; Rnext]
delta=(Rmax-1) * 0.5;       % maximum displacement of the initial positions. delta<=(Rmax-1)/2 preserves all the links
MaxSensingRadius=3;         % sensing radius of the agents (R_s)

%% Simulation parameters
Simulation=struct();
Simulation.Tmax =   20;     % maximum simulation time (simulation is stopped earlier if steady state is reached)
Simulation.deltaT = 0.25;   % sampling time step
Simulation.dT =     0.01;   % integration time step
Simulation.arena =  3;      % size of the simulation window
Simulation.drawON=false;    % draw swarm during simulation (if N is large slows down the simulation)
Simulation.drawTraj=true;   % draw trajectories of the agents (if N is large slows down the simulation)
Simulation.getMetrics=true; % acquire metrics during the simulation (getMetrics=false discard settling times and stop times)

%% Dynamic model of the agents
Dynamics=struct('model','FirstOrder', 'sigma',0, 'vMax', inf);
%Dynamics=struct('model','SecondOrder', 'sigma',0.1, 'vMax', inf);
%Dynamics=struct('model','CoupledSDEs', 'rateSpeed', 1, 'avgSpeed', avgSpeed0, 'sigmaSpeed', 1, 'rateOmega', 1, 'sigmaOmega', @(x)2*max(1-x/3,0), 'omega', zeros(N,1));
%Dynamics=struct('model','LevyWalk', 'alpha',0.005, 'sigma', 0.25);

%% Global and Local interaction functions
%GlobalIntFunction=struct('function','Lennard-Jones','parameters',[0.5, (D-1)*12], 'MaxSensingRadius', MaxSensingRadius, 'Gain', 1);
GlobalIntFunction=struct('function','PowerLaw-FiniteCutoff','parameters',[1, Rmax], 'MaxSensingRadius', MaxSensingRadius, 'Gain', 0.5);
%GlobalIntFunction=struct('function','Spears','parameters', [2 35]);  %from Spears2004
%GlobalIntFunction=struct('function','Morse','parameters',[0.2, 2]);
%GlobalIntFunction=struct('function','Modified-LJ','parameters',[]);  %from Torquato2009
%GlobalIntFunction=struct('function','None');

%LocalIntFunction=struct('function','Linear', 'LinkNumber',LinkNumber, 'DistanceRange', [0.6, 1.1], 'Gain', G_normal);
LocalIntFunction=struct('function','None', 'DistanceRange', [0, Rmax]);
%LocalIntFunction=struct('function','None');

smoothing = false;      % smooth temporal data with moving average



