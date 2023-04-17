function [xVec, vVec, stopTime] = Simulator(x0, v0, Tmax, Dynamics, drawON, getMetrics, GlobalIntFunction, LocalIntFunction)
%
%Simulator Executes a complete simulation of the swarm.
%   This function is called by a launcher script (Launcher, BruteForceTuning, ...).
%
%   [T_r, success, final_e_theta, final_e_L, final_e_d, finalGRadial, finalGNormal, stopTime] =...
%   ... Simulator(x0, LinkNumber, GainRadialDefault, GainNormalDefault, regularity_thresh, compactness_thresh,...
%   ... Tmax, sigma, drawON, getMetrics, IntFunctionStruct, AgentsRemoval, NoiseTest, MaxSensingRadius, alpha, beta, dynamicLattice)
%
%   Inputs:
%       x0 are the initial positions of the agents (Nx2 matrix)
%       LinkNumber is the desired number of links per agent (6=triangular
%           lattice, 4=square lattice, 3=hexagonal lattice) (scalar)
%       GainRadialDefault is the default value of G_radial (scalar)
%       GainNormalDefault is the default value of G_normal (scalar)
%       regularity_thresh is the threshold value for regularity metrics (e^*_theta) (scalar)
%       compactness_thresh is the threshold value for compactness metrics (e^*_L) (scalar)
%       Tmax is the maximum simulation time (scalar)
%       sigma is the standard deviation of noise (scalar)
%       drawON draw swarm during simulation (bool)
%       getMetrics acquire metrics during the simulation (getMetrics=false
%           discard settling times and stop times) (bool)
%       IntFunctionStruct description and parameters of the radial
%           interaction function (struct)
%       MaxSensingRadius is the sensing radius of the agents (R_s)
%       alpha and beta are the adaptation gains for the adaptive control (scalar)
%       dynamicLattice change lattice during the simulation (bool)
%       AgentsRemoval randomly remove agents during the simulation (bool)
%
%   Outputs:
%       T_r is the settling time (scalar)
%       success is true if the metrics are below the respective thresholda (bool)
%       final_e_theta final value of the regularity metrics (scalar)
%       final_e_L final value of the compactness metrics (scalar)
%       final_e_d final value of the lenght metrics (scalar)
%       finalGRadial final value of the radial gain (scalar)
%       finalGNormal final value of the normal gain (scalar)
%       stopTime final simulation instant (scalar)
%
%   See also: Launcher, BruteForceTuning, SequentialLauncher
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%

% %% Check input parameters
% assert(LinkNumber==6 | LinkNumber==4, "LinkNumber must be equal to 4 (square lattice) or 6 (triangular lattice)")
% assert(size(x0,2)==2, "x0 must be a Nx2 matrix")
% assert(regularity_thresh>0, "regularity_thresh must be a positive number")
% assert(compactness_thresh>0, "compactness_thresh must be a positive number")
% assert(Tmax>=0, "Tmax must be a non negative number")
% assert(sigma>=0, "sigma must be a non negative number")
% assert(MaxSensingRadius>=0, "MaxSensingRadius must be a non negative number")
% assert(islogical(drawON), "drawON must be a logical value")
% assert(islogical(getMetrics), "getMetrics must be a logical value")
% assert(islogical(dynamicLattice), "dynamicLattice must be a logical value")
% assert(islogical(AgentsRemoval), "AgentsRemoval must be a logical value")

%% Instantiate Simulation Window
Max = 10;   % amplitude of the simulation plane
Min = -Max;

if drawON
    figure;
    axis('equal',[Min Max Min Max])
    yticks([-10 -5 0 5 10])
    xticks([-10 -5 0 5 10])
    set(gca,'FontSize',14)
    hold on
end

vMax = 5; % maximum speed of the agents

SensingNumber = inf;    % max number of neighbours to interact with


%% simulation parameters
InteractionFactor = 1; % fraction of agents to interact with ]0,1] 
%(if InteractionFactor<1 a subset of agents is randomly selected by each agent at each update step to interact with)
if (InteractionFactor~=1); warning("InteractionFactor is NOT set to 1"); end

deltaT = 0.01;      % forward Euler integration step
deltaSample = 0.25; % time step for metrics acquisition
screenTimes=[0 Tmax/2-deltaT Tmax/2 Tmax]; % specify time instants to get simulation frames

% % steady state detection
% stopSamples=ceil(10/deltaSample);               % number of consecutive samples to detect steady state 
% regularity_ss_thresh=regularity_thresh/10;      % steady state detection threshold for the regularity metrics
% compactness_ss_thresh=compactness_thresh/10;    % steady state detection threshold for the compactness metrics

%% Inizialization
x=x0;
N=size(x,1);


%% Preallocate variables
count=0;                                    % sampling iteration
TSample = 0:deltaSample:Tmax;               % sampling time instants
xVec=nan([size(TSample,1)+1,size(x0)]);         % positions of the swarm
vVec=nan([size(TSample,1)+1,size(x0)]);         % velocity of the swarm

xVec(1,:,:)=x0;
vVec(1,:,:)=v0;
v=v0;

disp(['- Simulating ',GlobalIntFunction.function,' N=',num2str(N)])

%% Run Simulation
stopCondition=false;
t=0;
stopTime=nan;

while t<=Tmax && ~stopCondition

    
    % Compute Control Actions
    %[v, links, ~, G_radial, G_normal] = VFcontroller(x, G_radial, G_normal, min(RMax,MaxSensingRadius), RMin, SensingNumber, InteractionFactor, LinkNumber, deltaT, DeadZoneThresh, GlobalIntFunction, spin, MaxSensingRadius, alpha, beta);
    [vInteractions] = VFcontroller(x, SensingNumber, InteractionFactor, deltaT, GlobalIntFunction, LocalIntFunction);
    
    % Simulate Agents' Dynamics
    %x = SingleIntegrator(x, v, deltaT, vMax);
    [x, v] = Integrate(x, vInteractions, Dynamics, deltaT, vMax);
    
    if t>=TSample(count+1)
        count= count+1;
        
        xVec(count+1,:,:)=x;
        vVec(count+1,:,:)=v;
            
        if getMetrics || t>Tmax-2
            
            
%             % steady-state detection (with vibration exclusion)
%             if(count>stopSamples)
%                 if (abs(e_L(count-stopSamples:count,1)-e_L(count,1))< compactness_ss_thresh...
%                         | abs(e_L(count-stopSamples:count,1)-e_L(count-1,1))< compactness_ss_thresh)...
%                         & (abs(e_theta(count-stopSamples:count,1)-e_theta(count,1))< regularity_ss_thresh...
%                         | abs(e_theta(count-stopSamples:count,1)-e_theta(count-1,1))< regularity_ss_thresh)...
%                         & abs(G_normal_vec(count-stopSamples:count,1)-G_normal_vec(count,1))< G_normal_vec(count,1)*0.03
%                     stopCondition=true;
%                     stopTime=t;
%                 end
%             end
            
        end
        
        % plot swarm
        if drawON
            plotTrajectory(xVec, false, 'b');
            if isfield(LocalIntFunction, 'DistanceRange')
                plotSwarm(x, [], t, LocalIntFunction.DistanceRange(1), LocalIntFunction.DistanceRange(2), true, ones(N,1));
            else
                plotSwarm(x, [], t, inf, inf, true, ones(N,1));
            end
        end

    end
    
    
    t=t+deltaT;
end

%% PLOTS

% plot swarm
if drawON
    plotTrajectory(xVec, false, 'b');
    if isfield(LocalIntFunction, 'DistanceRange')
        plotSwarm(x, [], t, LocalIntFunction.DistanceRange(1), LocalIntFunction.DistanceRange(2), false, ones(N,1));
    else
        plotSwarm(x, [], t, inf, inf, false, ones(N,1));
    end
end



end

