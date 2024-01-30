function [xVec, uVec, vVec] = Simulator(x0, v0, Simulation, Dynamics, GlobalIntFunction, LocalIntFunction, Environment)
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
    Environment         struct = struct()
end

assert(ismember(size(x0,2), [2,3]), 'x0 must have second dimension equal to 2 or 3!')
assert(all(size(v0,2)==size(x0,2)), 'v0 must have same dimensions of x0!')

assert(Simulation.recordVideo <= Simulation.drawON, 'Simulation.recordVideo must be false if Simulation.drawON is false')

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

if length(Simulation.arena)==1
    Simulation.arena = [Simulation.arena, Simulation.arena];
end
    
%% Instantiate Simulation Window

if Simulation.drawON
    figure
    if isfield(LocalIntFunction, 'DistanceRange')
        plotSwarmInit(x0, 0, LocalIntFunction.DistanceRange(1), LocalIntFunction.DistanceRange(2), Simulation.arena, 1, false, false, true);
    else
        plotSwarmInit(x0, 0, inf, inf, Simulation.arena, 1, false, false, true);
    end
    
    if isfield(Environment,'Inputs') && isfield(Environment.Inputs,'Points')
        x_vec = linspace(-Simulation.arena(1)/2,Simulation.arena(1)/2,100);
        y_vec = linspace(-Simulation.arena(2)/2,Simulation.arena(2)/2,100);
        [x_mesh, y_mesh] = meshgrid(x_vec, y_vec);
        %F=scatteredInterpolant(Environment.Inputs.Points, Environment.Inputs.Values, 'linear', 'nearest');
        F = griddedInterpolant(Environment.Inputs.Points,Environment.Inputs.Values, 'linear', 'nearest');
        cmap = [linspace(1,1,100)',linspace(1,0.5,100)',linspace(1,0.5,100)'];
    end
end

%% Inizialization
epsilon = Simulation.dT/100;
x=x0;
N=size(x,1);


%% Preallocate variables
TSample = [0:Simulation.deltaT:Simulation.Tmax]';   % sampling time instants
xVec=nan([size(TSample,1),size(x0)]);               % positions of the swarm
vVec=nan([size(TSample,1),size(x0)]);               % velocities of the swarm
uVec=nan([size(TSample,1),size(x0)]);               % inputs of the swarm
envInput = zeros(N,1);

% xVec(1,:,:)=x0;
v=v0;

if Simulation.recordVideo
    video = VideoWriter('./Output/video','MPEG-4');
    video.FrameRate = 1/Simulation.deltaT/4;
    open(video);
%     currFrame = getframe(gcf);
%     writeVideo(video,currFrame);
end

disp(['- Simulating ',num2str(N),' ',Dynamics.model, ' agents in ', num2str(size(x0,2)),'D space with ',GlobalIntFunction.function,' interaction'])

%% Run Simulation
t=0;
count=1;                                        % sampling iteration

while t<Simulation.Tmax
    % Compute Control Actions from interactions
    forces = VFcontroller(x, GlobalIntFunction, LocalIntFunction, Simulation.dT, Simulation.InteractionFactor);
    
    % Compute environmental input
    if isfield(Environment,'Inputs')
        if isfield(Environment.Inputs,'Points')
            %envInput = interpn(Environment.Inputs.Points{1}, Environment.Inputs.Points{2}, Environment.Inputs.Values, x(:,1),x(:,2), linear, nearest);
            envInput = F(x(:,1),x(:,2));
        else
            envInput = ones(N,1) * interp1(Environment.Inputs.Times, Environment.Inputs.Values, t);
        end
    end
    
    % Acquire data
    if t>=TSample(count)-epsilon
        t=Simulation.deltaT*round(t/Simulation.deltaT);
        
        xVec(count,:,:)=x;
        vVec(count,:,:)=v;
        uVec(count,:,:)=forces;
        
        count= count+1;
        
        % plot swarm
        if Simulation.drawON
            cla
            hold on
            if isfield(Environment,'Inputs') && isfield(Environment.Inputs,'Points')
                imagesc(x_vec,y_vec,F(x_mesh',y_mesh')')
                colormap(cmap)
            end
            if Simulation.drawTraj; plotTrajectory(xVec, false, [0,0.7,0.9], Simulation.drawTraj); end
            if isfield(LocalIntFunction, 'DistanceRange')
                plotSwarm(x, [], t, LocalIntFunction.DistanceRange(1), LocalIntFunction.DistanceRange(2), true);
            else
                plotSwarm(x, [], t, inf, inf, true);
            end
            
            if Simulation.recordVideo
                currFrame = getframe(gcf);
                writeVideo(video,currFrame);
            end
        end
    end
    
    % Integrate Agents' Dynamics
    [x, v, Dynamics] = integrateAgents(x, v, forces, Dynamics, Simulation.dT, envInput);
    
    % Update time
    t=t+Simulation.dT;
    t=Simulation.dT*round(t/Simulation.dT);
end

xVec(count,:,:)=x;
vVec(count,:,:)=v;
uVec(count,:,:)=forces;

if Simulation.recordVideo
    close(video);
end

%% PLOTS

% plot swarm
if Simulation.drawON
    cla
    if isfield(Environment,'Inputs') && isfield(Environment.Inputs,'Points')
        imagesc(x_vec,y_vec,F(x_mesh',y_mesh')')
        colormap(cmap)
    end
    if Simulation.drawTraj; plotTrajectory(xVec, false, [0,0.7,0.9], Simulation.drawTraj); end
    if isfield(LocalIntFunction, 'DistanceRange')
        plotSwarm(squeeze(xVec(end,:,:)), [], t, LocalIntFunction.DistanceRange(1), LocalIntFunction.DistanceRange(2), false);
    else
        plotSwarm(squeeze(xVec(end,:,:)), [], t, inf, inf, false);
    end
end

end

