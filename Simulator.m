function [xVec, vVec, stopTime] = Simulator(x0, v0, Simulation, Dynamics, GlobalIntFunction, LocalIntFunction)
%
%Simulator Executes a complete simulation of the swarm.
%   This function is called by a launcher script (Launcher, BruteForceTuning, ...).
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

if Simulation.drawON
    figure;
    axis('equal',[Min Max Min Max])
    yticks([-10 -5 0 5 10])
    xticks([-10 -5 0 5 10])
    set(gca,'FontSize',14)
    hold on
end

SensingNumber = inf;    % max number of neighbours to interact with


%% simulation parameters
InteractionFactor = 1; % fraction of agents to interact with ]0,1] 
%(if InteractionFactor<1 a subset of agents is randomly selected by each agent at each update step to interact with)
if (InteractionFactor~=1); warning("InteractionFactor is NOT set to 1"); end

% % steady state detection
% stopSamples=ceil(10/deltaSample);               % number of consecutive samples to detect steady state 
% regularity_ss_thresh=regularity_thresh/10;      % steady state detection threshold for the regularity metrics
% compactness_ss_thresh=compactness_thresh/10;    % steady state detection threshold for the compactness metrics

%% Inizialization
x=x0;
N=size(x,1);


%% Preallocate variables
count=0;                                    % sampling iteration
TSample = 0:Simulation.deltaT:Simulation.Tmax;               % sampling time instants
xVec=nan([size(TSample,1)+1,size(x0)]);         % positions of the swarm
vVec=nan([size(TSample,1)+1,size(x0)]);         % velocity of the swarm

xVec(1,:,:)=x0;
vVec(1,:,:)=v0;
v=v0;

disp(['- Simulating ',Dynamics.model, ' with ',GlobalIntFunction.function,' interaction, N=',num2str(N)])

%% Run Simulation
stopCondition=false;
t=0;
stopTime=nan;

while t<=Simulation.Tmax && ~stopCondition

    
    % Compute Control Actions
    %[v, links, ~, G_radial, G_normal] = VFcontroller(x, G_radial, G_normal, min(RMax,MaxSensingRadius), RMin, SensingNumber, InteractionFactor, LinkNumber, deltaT, DeadZoneThresh, GlobalIntFunction, spin, MaxSensingRadius, alpha, beta);
    forces = VFcontroller(x, SensingNumber, InteractionFactor, Simulation.dT, GlobalIntFunction, LocalIntFunction);
    
    % Simulate Agents' Dynamics
    %x = SingleIntegrator(x, v, deltaT, vMax);
    [x, v, Dynamics] = Integrate(x, v, forces, Dynamics, Simulation.dT);
    
    if t>=TSample(count+1)
        count= count+1;
        
        xVec(count+1,:,:)=x;
        vVec(count+1,:,:)=v;
            
        if Simulation.getMetrics || t>Tmax-2
            
            
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
        if Simulation.drawON
            if Simulation.drawTraj; plotTrajectory(xVec, false, [0,0.7,0.9]); end
            if isfield(LocalIntFunction, 'DistanceRange')
                plotSwarm(x, [], t, LocalIntFunction.DistanceRange(1), LocalIntFunction.DistanceRange(2), true, ones(N,1));
            else
                plotSwarm(x, [], t, inf, inf, true, ones(N,1));
            end
        end

    end
    
    
    t=t+Simulation.dT;
end

%% PLOTS

% plot swarm
if Simulation.drawON
    if Simulation.drawTraj; plotTrajectory(xVec, false, [0,0.7,0.9]); end
    if isfield(LocalIntFunction, 'DistanceRange')
        plotSwarm(x, [], t, LocalIntFunction.DistanceRange(1), LocalIntFunction.DistanceRange(2), false, ones(N,1));
    else
        plotSwarm(x, [], t, inf, inf, false, ones(N,1));
    end
end



end

