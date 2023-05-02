%
%defaultParam Set the default values of the parameters.
%
%   See also: Launcher, SequentialLauncher
%   
%   Authors:    Andrea Giusti
%   Date:       2023
%

%% Default Parameters

N=50;          %number of agents (N)
LinkNumber=4;   %number of links (6=triangular lattice, 4=square lattice, 3=hexagonal lattice) (L)

% control gains
G_radial= 1;   % default value for square lattice 15 (G_r)
G_normal = 0;   % default value for square lattice  8 (G_n)

MaxSensingRadius=inf;   % sensing radius of the agents (R_s)

Simulation=struct();
Simulation.Tmax =   10;     % maximum simulation time (simulation is stopped earlier if steady state is reached)
Simulation.deltaT = 0.5;    % sampling time step
Simulation.dT =     0.01;   % integration time step
Simulation.drawON=true;    % draw swarm during simulation (if N is large slows down the simulation)
Simulation.drawTraj=true;  % draw trajectories of the agents (if N is large slows down the simulation)
Simulation.getMetrics=true; % acquire metrics during the simulation (getMetrics=false discard settling times and stop times)


% descrption of the radial and normal interaction functions
GlobalIntFunction=struct('function','Lennard-Jones','parameters',[0.5, 12], 'MaxSensingRadius', MaxSensingRadius, 'Gain', G_radial);
%GlobalIntFunction=struct('function','Spears','parameters', [2 35]);  %from Spears2004
%GlobalIntFunction=struct('function','Morse','parameters',[0.2, 2]);
%GlobalIntFunction=struct('function','Modified-LJ','parameters',[]);  %from Torquato2009
%GlobalIntFunction=struct('function','None');

%LocalIntFunction=struct('function','Linear', 'LinkNumber',LinkNumber, 'DistanceRange', [0.6, 1.1], 'Gain', G_normal);
LocalIntFunction=struct('function','None', 'DistanceRange', [0, 1.1]);

% adaptation gains
alpha = 0;      
beta = 0;

% thresholds
regularity_thresh=0.2;      % threshold value for regularity metrics (e^*_theta)
compactness_thresh=0.3;     % threshold value for compactness metrics (e^*_L)

Rmax= (sqrt(3)+1)/2;    % maximum lenght of a link (R_a). Must be in [1; sqrt(3)]
delta=(Rmax-1) * 0.5;   % maximum displacement of the initial positions. delta<=(Rmax-1)/2 preserves all the links

% robustness tests
AgentsRemoval=false;
NoiseTest=false;
dynamicLattice = false;


