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

% Directory to save the results of the simulations.
% Set outputDir='' to prevent automatic saving.
%outputDir='./Output';
outputDir='';

N=10;                      %number of agents (N)
LinkNumber=6*(D-1);         %number of links per agent in the lattice configuration (L)
                            %If D=2 then 6=triangular lattice, 4=square lattice, 3=hexagonal lattice
                            %If D=3 then 6=cubic lattice, 12=thetradic-octaedric lattice
                            
smoothing = false;          % smooth temporal data with moving average

Rmax= (sqrt(5-D)+1)/2;      % maximum lenght of a link (R_a). Must be in [1; Rnext]
delta=(Rmax-1) * 0.5;       % maximum displacement of the initial positions. delta<=(Rmax-1)/2 preserves all the links
MaxSensingRadius=3;         % sensing radius of the agents (R_s)

%% Simulation parameters
% All these fields are mandatory
Simulation=struct();
Simulation.Tmax =   20;     % maximum simulation time (simulation is stopped earlier if steady state is reached)
Simulation.deltaT = 0.1;   % sampling time step
Simulation.dT =     0.01;   % integration time step
Simulation.arena =  800;    % size of the simulation window
Simulation.drawON=false;    % draw swarm during simulation (if N is large slows down the simulation)
Simulation.drawTraj=15;   % draw trajectories of the agents (if N is large slows down the simulation)
Simulation.recordVideo=true;% record video of the simulation (if true drawON must be true)
Simulation.getMetrics=true; % acquire metrics during the simulation (getMetrics=false discard settling times and stop times)

%% Dynamic model of the agents
% Initial velocities for CoupledSDEs and LevyWalk.
avgSpeed0   = 10;
sigmaSpeed0 = 0;

% These parameters are used in integrateAgents.

%Dynamics=struct('model','FirstOrder', 'sigma',0, 'vMax', inf);
%Dynamics=struct('model','SecondOrder', 'sigma',0.1, 'vMax', inf);
Dynamics=struct('model','IndependentSDEs', 'rateSpeed', 0.1, 'avgSpeed', avgSpeed0, 'sigmaSpeed', 3, 'rateOmega', 0.31, 'sigmaOmega', 0.23, 'omega', normrnd(0,0.3,N,1));
%Dynamics=struct('model','CoupledSDEs', 'rateSpeed', 1, 'avgSpeed', avgSpeed0, 'sigmaSpeed', 1, 'rateOmega', 1, 'sigmaOmega', @(x)2*max(1-x/3,0), 'omega', zeros(N,1));
%Dynamics=struct('model','LevyWalk', 'alpha',0.005, 'sigma', 0.15);

%% Global interaction function for long distance interactions
% These parameters are used in globalInteractionForce.

%GlobalIntFunction=struct('function','Lennard-Jones','parameters',[0.5, (D-1)*12], 'MaxSensingRadius', MaxSensingRadius, 'Gain', 1);
%GlobalIntFunction=struct('function','PowerLaw-FiniteCutoff','parameters',[1, Rmax], 'MaxSensingRadius', MaxSensingRadius, 'Gain', 0.5);
%GlobalIntFunction=struct('function','Spears','parameters', [2 35]);  %from Spears2004
%GlobalIntFunction=struct('function','Morse','parameters',[0.2, 2]);
%GlobalIntFunction=struct('function','Modified-LJ','parameters',[]);  %from Torquato2009
GlobalIntFunction=struct('function','None');

%% Local interaction function for short distance interactions
% These parameters are used in localInteractionForce.
% If LocalIntFunction has a 'DistanceRange' field it is used to compute and 
% plot the links between the agents.

%LocalIntFunction=struct('function','Linear', 'LinkNumber',LinkNumber, 'DistanceRange', [0.6, 1.1], 'Gain', 1);
%LocalIntFunction=struct('function','None', 'DistanceRange', [0, Rmax]);
LocalIntFunction=struct('function','None');

% Set an optional rotation matrix to apply non-radial local forces.
% Normal intercations can be used to form square lattices (only in 2D).
% LocalIntFunction.Rotation = [0 1; -1 0];  % 90deg rotation matrix (optional)


%% Add subfolders to the Matlab path
current_folder = fileparts(which('defaultParam'));
addpath(genpath(current_folder));

