%
%defaultParam Set the default values of the parameters.
%
%   See also: Launcher, SequentialLauncher
%   
%   Authors:    Andrea Giusti
%   Date:       2023
%

%% Default Parameters

N=100;          %number of agents (N)
LinkNumber=6;   %number of links (6=triangular lattice, 4=square lattice, 3=hexagonal lattice) (L)

% descrption of the radial interaction function
RadialIntFunction=struct('function','Lennard-Jones','parameters',[0.5, 12]);
%RadialIntFunction=struct('function','Spears','parameters', [2 35]);  %from Spears2004
%RadialIntFunction=struct('function','Morse','parameters',[0.2, 2]);
%RadialIntFunction=struct('function','Modified-LJ','parameters',[]);  %from Torquato2009

% control gains
G_radial= 1;   % default value for square lattice 15 (G_r)
G_normal = 0;   % default value for square lattice  8 (G_n)

% adaptation gains
alpha = 0;      
beta = 0;

% thresholds
regularity_thresh=0.2;      % threshold value for regularity metrics (e^*_theta)
compactness_thresh=0.3;     % threshold value for compactness metrics (e^*_L)

Tmax=10;    % maximum simulation time (simulation is stopped earlier if steady state is reached)

sigma = 0;  % standard deviation of noise in agents' dynamics

MaxSensingRadius=inf;   % sensing radius of the agents (R_s)

Rmax= (sqrt(3)+1)/2;    % maximum lenght of a link (R_a). Must be in [1; sqrt(3)]
delta=(Rmax-1) * 0.5;   % maximum displacement of the initial positions. delta<=(Rmax-1)/2 preserves all the links

% robustness tests
AgentsRemoval=false;
NoiseTest=false;
dynamicLattice = false;

%output options
drawON=false;       % draw swarm during simulation (if N is large slows down the simulation)
getMetrics=true;    % acquire metrics during the simulation (getMetrics=false discard settling times and stop times)




