%
%RobotariumSimulator allows to run experiments in the Robotarium and
%   simulations in the Robotarium simulator.
%
%   At the beginning of the experiment agents go to inital positions
%   generated by the Robotarium.
%   Then they attract each other using a consenus algorithm, so to form a compact group.
%   Then the lattice formation algorithm starts.
%
%   Notes: 
%       To run this script you need to download the Robotarium source code
%       form https://www.robotarium.gatech.edu/downloads and add it to your
%       Matlab path.
%       For more details and to use the Robotarium refer to https://www.robotarium.gatech.edu.
%
%   See also: Launcher
%
%   Author:         Paul Glotfelter
%   Date:           10/04/2016
%
%   Modified by:    Andrea Giusti and Gian Carlo Maffettone
%   Date:           2022
%

%% Clear environment
close all
clc 
clear

%% Parameters
N=20;           % number of agents, max 20 (N)
LinkNumber=4;   %number of links (6=triangular lattice, 4=square lattice) (L)

dynLattice = false; % change lattice during the experiment

% descrption of the radial interaction function
IntFunctionStruct=struct('function','Lennard-Jones','parameters',[0.15, 5]);

% control gains
G_radial= 0.8;
G_normal = 0.4;
G_radial_vec = G_radial* ones(N,1);
G_normal_vec = G_normal* ones(N,1);

% adaptation gains
alpha = 0;      
beta = 0;

% thresholds
regularity_thresh=0.2;              % threshold value for regularity metrics (e^*_theta)
compactness_thresh=0.3;             % threshold value for compactness metrics (e^*_L)

sigma = 0;                          % standard deviation of noise

MaxSensingRadius=inf;               % sensing radius of the agents (R_s)
RMax = 1.1;                         % maximum distance of adjacent agents
RMin = 0.6;                         % minimum distance of adjacent agents
SensingNumber = inf;                % max number of neighbours to interact with
DeadZoneThresh=regularity_thresh;   % amplitude of the dead zone in gain adaptation law
InteractionFactor = 1;              % fraction of agents to interact with [0,1] 
%(if InteractionFactor<1 a subset of agents is randomly selected by each agent at each update step to interact with)

deltaT = 0.033;     % integration time step

s=ones(N, 1);       %Spears' spin
if strcmp(IntFunctionStruct.function,'Spears') && (LinkNumber==4 || LinkNumber == 3)
    s(1:2:N, 1)=0;
end

ScaleFactor = 2;        % to map our environment in Robotarium's environment
K = 1;                  % aggregation gain
t_aggregation = 1000;   % iterations for initail aggregation

WheelVelocityLimit = 12.5/2;    %max velocity for unicycle default 12.5
Vmax = 0.1;                     %max velocity for SI dynamics (robotarium units)
Vmin = 0;                       %min velocity for SI dynamics (robotarium units)

SafetyRad=0.11;

% create initial conditions
rng(4,'twister'); % reproducible results
initial_positions = generate_initial_conditions(N, 'Spacing', 0.3, 'Width', 2.8, 'Height', 1.5);

% robotarium object
r = Robotarium('NumberOfRobots', N, 'ShowFigure', true, 'InitialConditions', initial_positions);

% set collision avoidance policy
si_to_uni_dynamics = create_si_to_uni_dynamics_with_backwards_motion('AngularVelocityLimit', pi/10);
uni_barrier_certificate = create_uni_barrier_certificate_with_boundary('SafetyRadius', SafetyRad, 'BarrierGain', 50, 'BoundaryPoints', [-1.4, 1.4, -0.8, 0.8], 'WheelVelocityLimit', WheelVelocityLimit, 'ProjectionDistance', 0.03);

% Select the number of iterations for the experiment.
iterations = 5000;

%setup for data saving
positions = zeros(N,2,iterations);
velocities = zeros(N,2,iterations);

% Plot the iteration and time in the lower left. Note when run on your 
% computer, this time is based on your computers simulation time. For a
% better approximation of experiment time on the Robotarium when running
% this simulation on your computer, multiply iteration by 0.033. 
% Marker, font, and line sizes
marker_size_goal = determine_marker_size(r, 0.20);
marker_size_robot = determine_robot_marker_size(r);
font_size = determine_font_size(r, 0.05);
line_width = 5;
start_time = tic; %The start time to compute time elapsed.
iteration_caption = sprintf('Iteration %d', 0);
time_caption = sprintf('Total Time Elapsed %0.2f', toc(start_time));
iteration_label = text(-1.5, -0.8, iteration_caption, 'FontSize', font_size, 'Color', 'r', 'FontWeight', 'bold');
time_label = text(-1.5, -0.9, time_caption, 'FontSize', font_size, 'Color', 'r', 'FontWeight', 'bold');


% Iterate for the previously specified number of iterations
for t = 1:iterations
    
    if dynLattice && t > (iterations - t_aggregation)/3 + t_aggregation && t <= 2*(iterations - t_aggregation)/3 + t_aggregation
        LinkNumber = 6;
    else
        LinkNumber = 4;
    end
    
    
    % Retrieve the most recent poses from the Robotarium.  The time delay is
    % approximately 0.033 seconds
    x = r.get_poses();
    
    %import data and convert from robotarium units
    positions(:,:,t) = x(1:2,:)' * ScaleFactor;
    
    %% Compute control actions
        
    if t < t_aggregation        % aggregation
        dxi = -K * x(1:2,:)';
        
    else                        % lattice formation
        lines = drawLines(x(1:2,:)',RMin/ScaleFactor,RMax/ScaleFactor);
        
        % Compute Control Actions
        [dxi, ~, ~, ~, ~]= VFcontroller(x(1:2,:)'*ScaleFactor, G_radial_vec, G_normal_vec, min(RMax,MaxSensingRadius), RMin, SensingNumber, InteractionFactor, LinkNumber, deltaT, DeadZoneThresh, IntFunctionStruct, s, MaxSensingRadius, alpha, beta);
    end
    
    dxi = dxi'/ScaleFactor; % conversion to robotarium units

    
    %% Send velocities to agents
    
    %dxi = si_barrier_certificate(dxi, x(1:2, :));   %apply control barrier function to single integrator dynamics
    dxu = si_to_uni_dynamics(dxi, x);                %convert control action to unicycle commands
    dxu = uni_barrier_certificate(dxu, x);           %apply control barrier function to unicycle dynamics (prevents actuator limits errors)

    velocities(:,:,t) = dxu';
    
    r.set_velocities(1:N, dxu);
        
    r.step();
    
    if t >= t_aggregation
        delete(lines)
    end
    
    %% Update Iteration and Time marker
    iteration_caption = sprintf('Iteration %d', t);
    time_caption = sprintf('Total Time Elapsed %0.2f', toc(start_time));
    
    iteration_label.String = iteration_caption;
    time_label.String = time_caption;
end

lines = drawLines(x(1:2,:)',RMin/ScaleFactor,RMax/ScaleFactor);


% datasaving
simulation = struct('positions', positions, 'velocities', velocities, 'LinkNumber', LinkNumber, 'dynLattice', dynLattice);
save('simulation.mat', 'simulation');


% We should call r.call_at_scripts_end() after our experiment is over!
r.debug();


%% Helper Functions

% Marker Size Helper Function to scale size of markers for robots with figure window
% Input: robotarium class instance
function marker_size = determine_robot_marker_size(robotarium_instance)

% Get the size of the robotarium figure window in pixels
curunits = get(robotarium_instance.figure_handle, 'Units');
set(robotarium_instance.figure_handle, 'Units', 'Pixels');
cursize = get(robotarium_instance.figure_handle, 'Position');
set(robotarium_instance.figure_handle, 'Units', curunits);

% Determine the ratio of the robot size to the x-axis (the axis are
% normalized so you could do this with y and figure height as well).
robot_ratio = (robotarium_instance.robot_diameter + 0.03)/...
    (robotarium_instance.boundaries(2) - robotarium_instance.boundaries(1));

% Determine the marker size in points so it fits the window. cursize(3) is
% the width of the figure window in pixels. (the axis are
% normalized so you could do this with y and figure height as well).
marker_size = cursize(3) * robot_ratio;

end

% Marker Size Helper Function to scale size with figure window
% Input: robotarium instance, desired size of the marker in meters
function marker_size = determine_marker_size(robotarium_instance, marker_size_meters)

% Get the size of the robotarium figure window in pixels
curunits = get(robotarium_instance.figure_handle, 'Units');
set(robotarium_instance.figure_handle, 'Units', 'Pixels');
cursize = get(robotarium_instance.figure_handle, 'Position');
set(robotarium_instance.figure_handle, 'Units', curunits);

% Determine the ratio of the robot size to the x-axis (the axis are
% normalized so you could do this with y and figure height as well).
marker_ratio = (marker_size_meters)/(robotarium_instance.boundaries(2) -...
    robotarium_instance.boundaries(1));

% Determine the marker size in points so it fits the window. cursize(3) is
% the width of the figure window in pixels. (the axis are
% normalized so you could do this with y and figure height as well).
marker_size = cursize(3) * marker_ratio;

end

% Font Size Helper Function to scale size with figure window
% Input: robotarium instance, desired height of the font in meters
function font_size = determine_font_size(robotarium_instance, font_height_meters)

% Get the size of the robotarium figure window in point units
curunits = get(robotarium_instance.figure_handle, 'Units');
set(robotarium_instance.figure_handle, 'Units', 'Pixels');
cursize = get(robotarium_instance.figure_handle, 'Position');
set(robotarium_instance.figure_handle, 'Units', curunits);

% Determine the ratio of the font height to the y-axis
font_ratio = (font_height_meters)/(robotarium_instance.boundaries(4) -...
    robotarium_instance.boundaries(3));

% Determine the font size in points so it fits the window. cursize(4) is
% the hight of the figure window in points.
font_size = cursize(4) * font_ratio;

end
