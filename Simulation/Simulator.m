function [xVec, uVec] = Simulator(x0, v0, Simulation, Dynamics, GlobalIntFunction, LocalIntFunction)
%
%Simulator executes a complete simulation of the swarm.
%   This function is called by a launcher script (Launcher, SequentialLauncher...).
%
%   [xVec, uVec] = Simulator(x0, v0, Simulation, Dynamics, GlobalIntFunction, LocalIntFunction)
%
%   Inputs:
%       x0                  Initial positions of the agents     (NxD matrix)
%       v0                  Initial velocities of the agents    (NxD matrix)
%       Simulation          Simulation parameters               (struct)
%       Dynamics            Dynamics of the agents              (struct)
%       GlobalIntFunction   Long distance interaction           (struct = struct('function','None'))
%       LocalIntFunction    Shorts distance interaction         (struct = struct('function','None'))
%
%   Outputs:
%       xVec                Positions of the agents             (TIMExNxD matrix)
%       uVec                Virtual forces acting on the agents (TIMExNxD matrix)
%
%   See also: Simulator, Launcher
%
%   Authors:    Andrea Giusti
%   Date:       2023
%

%% Validate input arguments
arguments
    x0                  double
    v0                  double
    Simulation          struct
    Dynamics            struct
    GlobalIntFunction   struct = struct('function','None')
    LocalIntFunction    struct = struct('function','None')
end

assert(ismember(size(x0,2), [2,3]), 'x0 must have second dimension equal to 2 or 3!')
assert(all(size(v0,2)==size(x0,2)), 'v0 must have same dimensions of x0!')

if ~isfield(Simulation, 'InteractionFactor')
    Simulation.InteractionFactor = 1; % fraction of agents to interact with ]0,1] 
end
%(if InteractionFactor<1 a subset of agents is randomly selected by each agent at each update step to interact with)
assert(Simulation.InteractionFactor <=1 & Simulation.InteractionFactor >0, 'Simulation.InteractionFactor must be in range ]0;1]')
if (Simulation.InteractionFactor~=1); warning("Simulation.InteractionFactor is NOT set to 1"); end

if ~isfield(GlobalIntFunction, 'SensingNumber')
    GlobalIntFunction.SensingNumber = inf; % fraction of agents to interact with ]0,1] 
end
assert(GlobalIntFunction.SensingNumber>0, 'GlobalIntFunction.SensingNumber must be positive')
if (GlobalIntFunction.SensingNumber~=inf); warning("SensingNumber is NOT set to inf"); end

%% Instantiate Simulation Window

if Simulation.drawON
    if isfield(LocalIntFunction, 'DistanceRange')
        plotSwarmInit(x0, 0, LocalIntFunction.DistanceRange(1), LocalIntFunction.DistanceRange(2), Simulation.arena, 1, false, false, true);
    else
        plotSwarmInit(x0, 0, inf, inf, Simulation.arena, 1, false, false, true);
    end
end

%% Inizialization
x=x0;
N=size(x,1);


%% Preallocate variables
count=0;                                        % sampling iteration
TSample = 0:Simulation.deltaT:Simulation.Tmax;  % sampling time instants
xVec=nan([size(TSample,1)+1,size(x0)]);         % positions of the swarm
uVec=nan([size(TSample,1)+1,size(x0)]);         % inputs of the swarm

xVec(1,:,:)=x0;
v=v0;

disp(['- Simulating ',num2str(N),' ',Dynamics.model, ' agents in ', num2str(size(x0,2)),'D space with ',GlobalIntFunction.function,' interaction'])

%% Run Simulation
t=0;

while t<=Simulation.Tmax
    
    % Compute Control Actions
    forces = VFcontroller(x, GlobalIntFunction, LocalIntFunction, Simulation.dT, Simulation.InteractionFactor);
    
    % Simulate Agents' Dynamics
    [x, v, Dynamics] = integrateAgents(x, v, forces, Dynamics, Simulation.dT);
    
    if t>=TSample(count+1)
        count= count+1;
        
        xVec(count+1,:,:)=x;
        uVec(count+1,:,:)=forces;
        
        % plot swarm
        if Simulation.drawON
            if Simulation.drawTraj; plotTrajectory(xVec, false, [0,0.7,0.9]); end
            if isfield(LocalIntFunction, 'DistanceRange')
                plotSwarm(x, [], t, LocalIntFunction.DistanceRange(1), LocalIntFunction.DistanceRange(2), true);
            else
                plotSwarm(x, [], t, inf, inf, true);
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
        plotSwarm(x, [], t, LocalIntFunction.DistanceRange(1), LocalIntFunction.DistanceRange(2), false);
    else
        plotSwarm(x, [], t, inf, inf, false);
    end
end

end

