function [xVec, uVec, vVec] = Simulator(x0, v0, Simulation, Dynamics, Render, GlobalIntFunction, LocalIntFunction, Environment)
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
    Render              struct = struct()
    GlobalIntFunction   struct = struct('function','None')
    LocalIntFunction    struct = struct('function','None')
    Environment         struct = struct()
end

assert(ismember(size(x0,2), [2,3]), 'x0 must have second dimension equal to 2 or 3!')
assert(all(size(v0,2)==size(x0,2)), 'v0 must have same dimensions of x0!')

assert(Render.recordVideo <= Render.drawON, 'Render.recordVideo must be false if Render.drawON is false')

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

cmap = Render.cmap_inputs;
if isfield(Environment,'Inputs') && isfield(Environment.Inputs,'Points')
    x_vec = linspace(-Simulation.arena(1)/2,Simulation.arena(1)/2,1000);
    y_vec = linspace(-Simulation.arena(2)/2,Simulation.arena(2)/2,1000);
    [x_mesh, y_mesh] = meshgrid(x_vec, y_vec);
    
    if Environment.Inputs.Times(1) > 0
        Environment.Inputs.Points = {[0,1], [0,1]; Environment.Inputs.Points{:,:}};
        Environment.Inputs.Values = {zeros(2); Environment.Inputs.Values{:}};
        Environment.Inputs.Times = [0, Environment.Inputs.Times];
    end
end

%% Instantiate Simulation Window

if Render.drawON
    figure
    if isfield(Environment,'Inputs') && isfield(Environment.Inputs,'Points')
        F = griddedInterpolant(Environment.Inputs.Points(1,:),Environment.Inputs.Values{1}, 'linear', 'nearest');
        imagesc(x_vec,y_vec,F(x_mesh',y_mesh')')
        colormap(cmap)
    end
    if isfield(LocalIntFunction, 'DistanceRange')
        plotSwarmInit(x0, 0, LocalIntFunction.DistanceRange(1), LocalIntFunction.DistanceRange(2), Render.window, [Render.window(2)-Render.window(1), Render.window(4)-Render.window(3)]/2, false, false, true, Render.agentShape, Render.agentSize, Render.agentsColor, x0-v0);
    else
        plotSwarmInit(x0, 0, inf, inf, Render.window, [Render.window(2)-Render.window(1), Render.window(4)-Render.window(3)]/2, false, false, true, Render.agentShape, Render.agentSize, Render.agentsColor, x0-v0);
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
forces=zeros(size(x0));
envInput = zeros(N,1);

% xVec(1,:,:)=x0;
v=v0;

if Render.recordVideo
    SSdir = getSSfolder();
    out_dir = fullfile(SSdir,'Output');
    video = VideoWriter(fullfile(out_dir,'video'),'MPEG-4');
    video.FrameRate = Render.frameRate;
    open(video);
end

log_txt = ['- Simulating ',num2str(N),' ',Dynamics.model, ' agents in ', num2str(size(x0,2)),'D space'];
if ~strcmp(GlobalIntFunction.function, 'None'); log_txt = [log_txt,' with ',GlobalIntFunction.function,' interaction']; end
if isfield(Environment,'Inputs')
    if isfield(Environment.Inputs,'Points') && isfield(Environment.Inputs,'Times')
        log_txt = [log_txt,' with spatiotemporal inputs'];
    elseif isfield(Environment.Inputs,'Points')
        log_txt = [log_txt,' with spatial inputs'];
    else
        log_txt = [log_txt,' with temporal inputs'];
    end
end
disp(log_txt);

%% Run Simulation
t=0;
count=1;                                        % sampling iteration

while t<Simulation.Tmax
    % Compute inputs from interactions
    if ~strcmp(GlobalIntFunction.function, 'None') || ~strcmp(LocalIntFunction.function, 'None')
        forces = VFcontroller(x, GlobalIntFunction, LocalIntFunction, Simulation.dT, Simulation.InteractionFactor);
    end
    
    % Compute environmental inputs
    if isfield(Environment,'Inputs')
        if isfield(Environment.Inputs,'Points') && ~all(strcmp(Environment.Inputs.Points,'None'),'all')
            inputIndex = find(Environment.Inputs.Times <= t, 1, 'last');
            F = griddedInterpolant(Environment.Inputs.Points(inputIndex,:),Environment.Inputs.Values{inputIndex}, 'linear', 'nearest');
            envInput = F(x(:,1),x(:,2));
        elseif isfield(Environment.Inputs,'Times') && ~all(strcmp(Environment.Inputs.Times,'None'),'all')
            envInput = ones(N,1) * interp1(Environment.Inputs.Times, Environment.Inputs.Values, t, Environment.Inputs.InterpMethod);
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
        if Render.drawON
            cla
            hold on
            if isfield(Environment,'Inputs')
                if isfield(Environment.Inputs,'Points') && ~all(strcmp(Environment.Inputs.Points,'None'),'all')
                    imagesc(x_vec,y_vec,F(x_mesh',y_mesh')')
                elseif isfield(Environment.Inputs,'Times') && ~all(strcmp(Environment.Inputs.Times,'None'),'all')
                    set(gca,'Color',interp1(linspace(0,1,100),cmap,interp1(Environment.Inputs.Times, Environment.Inputs.Values, t, Environment.Inputs.InterpMethod)))
                end
            end
            if Render.drawTraj; plotTrajectory(xVec, false, [0,0.7,0.9], Render.drawTraj); end
            if isfield(Environment,'boundary'); plotBoundary(Environment.boundary); end
            if isfield(LocalIntFunction, 'DistanceRange')
                plotSwarm(x, t, LocalIntFunction.DistanceRange(1), LocalIntFunction.DistanceRange(2), ~Render.recordVideo, ones(size(x,1), 1), false, Render.agentShape, Render.agentSize, Render.agentsColor, x-v);
            else
                plotSwarm(x, t, inf, inf, ~Render.recordVideo, ones(size(x,1), 1), false, Render.agentShape, Render.agentSize, Render.agentsColor, x-v);
            end
            
            if Render.recordVideo
                currFrame = getframe(gcf);
                writeVideo(video,currFrame);
            end
        end
    end
    
    % Integrate Agents' Dynamics
    [x, v, Dynamics] = integrateAgents(x, v, forces, Dynamics, Simulation.dT, envInput);
    
    % Compute boundaries interaction
    if isfield(Environment,'boundary')
        [x, v, out_agents] = boundaryInteraction(x, v, Environment.boundary);
    end
    
    % Update time
    t=t+Simulation.dT;
    t=Simulation.dT*round(t/Simulation.dT);
end

xVec(count,:,:)=x;
vVec(count,:,:)=v;
uVec(count,:,:)=forces;



%% PLOTS

% plot swarm
if Render.drawON
    cla
    if isfield(Environment,'Inputs') && isfield(Environment.Inputs,'Points')
        imagesc(x_vec,y_vec,F(x_mesh',y_mesh')')
        %colormap(cmap)
    end
    if Render.drawTraj; plotTrajectory(xVec, false, [0,0.7,0.9], Render.drawTraj); end
    if isfield(LocalIntFunction, 'DistanceRange')
        plotSwarm(squeeze(xVec(end,:,:)), t, LocalIntFunction.DistanceRange(1), LocalIntFunction.DistanceRange(2), false, ones(size(x,1), 1), false, Render.agentShape, Render.agentSize, Render.agentsColor, x-v);
    else
        plotSwarm(squeeze(xVec(end,:,:)), t, inf, inf, false, ones(size(x,1), 1), false, Render.agentShape, Render.agentSize, Render.agentsColor, x-v);
    end
    if isfield(Environment,'boundary'); plotBoundary(Environment.boundary); end
    
    if Render.recordVideo
        currFrame = getframe(gcf);
        writeVideo(video,currFrame);
        close(video);
    end
end

end

