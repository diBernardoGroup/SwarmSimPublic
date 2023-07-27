function [v, links, angErr, GainRadial, GainNormal] = VFcontroller(x, oldGainRadial, oldGainNormal, RMax, RMin, SensingNumber, InteractionFactor, LinkNumber, deltaT, DeadZoneThresh, IntFunctionStruct, spin, MaxSensingRadius, alpha, beta, sigma_measure, compassBias)
%
%VFcontroller implements the virtual forces controller and computes the
%   velocities (control input) of the agents.
%   This function is called by Simulator.
%
%   [v, links, angErr, GainRadial, GainNormal] = ...
%   ... VFcontroller(x, oldGainRadial, oldGainNormal, RMax, RMin,...
%   ... SensingNumber, InteractionFactor, LinkNumber, deltaT, DeadZoneThresh,...
%   ... IntFunctionStruct, spin, MaxSensingRadius, alpha, beta) 
%
%   Inputs:
%       x are the positions of the agents (Nx2 matrix)
%       LinkNumber is the desired number of links per agent (scalar)
%       oldGainRadial is the last value of G_radial (vector)
%       oldGainNormal is the last value of G_normal (vector)
%       oldGainRand is the last value of G_random (vector)
%       RMax and RMin are the distances that define the adjacency set (scalar)
%       SensingNumber is the maximum number of neighbours agents can sense (integer)
%       InteractionFactor is the fraction of agents to interact with (sclar)
%       IntFunctionStruct description and parameters of the radial
%           interaction function (struct)
%       MaxSensingRadius is the sensing radius of the agents (scalar)
%       deltaT integration time step (scalar)
%       spin is the spin of the agents for Spears' algorithm (vector)
%       DeadZoneThresh amplitute of the adaptation deadhzone (scalar)
%       alpha and beta are the adaptation gains for the adaptive control (scalar)
%
%   Outputs:
%       v are the velocities of the agents (Nx2 matrix)
%       links are the number of links of each agent (vector)
%       angErr average angular error of each agent (vector)
%       GainRadial updated values of the radial gain (vector)
%       GainNormal updated values of the normal gain (vector)
%
%   See also: Simulator
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%

%% Check input parameters
assert(LinkNumber==6 | LinkNumber==4, "LinkNumber must be equal to 4 (square lattice) or 6 (triangular lattice)")
assert(size(x,2)==2, "x must be a Nx2 matrix")
assert(deltaT>=0, "deltaT must be a non negative number")
assert(RMax>=RMin, "RMax must be >=RMin")
assert(RMin>=0, "RMin must be a non negative number")
assert(ceil(SensingNumber)==floor(SensingNumber), "SensingNumber must be an integer number")
assert(InteractionFactor>0 && InteractionFactor<=1, "InteractionFactor must be in ]0; 1]")
assert(MaxSensingRadius>=0, "MaxSensingRadius must be a non negative number")
assert(DeadZoneThresh>=0, "DeadZoneThresh must be a non negative number")
assert(all(spin==0 | spin==1), "spin must be equal to 0 or 1")

%% Instantiate Variables
    N=size(x,1);
    v= zeros(N,2);
    links = zeros(N,1);
    angErr = zeros(N,1);
    GainRadial = zeros(N,1);
    GainNormal = zeros(N,1);
    R= [0 1; -1 0];         % 90deg rotation matrix
    
%% for each agent...
    for i=1:N
        %% select agents to interact with
        if SensingNumber < size(x,1)
            xOthers = GetNcloser(x(i,:), x, SensingNumber);
        else
            xOthers = x;
        end
        
        % randomly select a subset of agents to interact with (if InteractionFactor < 1)
        if InteractionFactor < 1
            index = randsample(size(xOthers,1),ceil(N*InteractionFactor));
            xRand = xOthers(index,:);
        else
            xRand = xOthers;
        end
        
        % compute the adjacency set
        [xNeighbours, NeigIndices] = getNeighbours(x(i,:), xRand, RMin, RMax);
        
        
        %% links and angular error
        links(i)=size(xNeighbours,1);
        angErr(i)=getAngleError(x(i,:), 0, xNeighbours, LinkNumber);
        
        %% compute gains with adaptation law
        [GainRadial(i), GainNormal(i)]=adaptationLaw(alpha, beta, oldGainRadial(i), oldGainNormal(i), LinkNumber-links(i), angErr(i), deltaT, DeadZoneThresh, LinkNumber);
        
        
        %% compute interactions with neighbours
        % compute the distances from the neighbours (||r_ij||)
        distances = vecnorm(x(i,:)-xRand, 2, 2); 
        distances = distances+randn(size(distances))*sigma_measure; % add measurement noise
        
        % Renormalize distances if using Spear'spin algorithm (sqaure or hexagonal lattice)
        if strcmp(IntFunctionStruct.function,'Spears') && (LinkNumber==4 || LinkNumber == 3)
            same_spin = find(spin == spin(i));
            if LinkNumber == 4
                distances(same_spin) = distances(same_spin)/sqrt(2);
            else
                distances(same_spin) = distances(same_spin)/sqrt(3);
            end
        end
        
        % compute angular errors (theta_ij^err)
        angErrNeigh = getAngularErrNeigh(x(i,:), 0, xNeighbours, LinkNumber);
        angErrNeigh = angErrNeigh+randn(size(angErrNeigh)).*sigma_measure*pi/LinkNumber;   % add measurement noise i.i.d. on each angular measure
        %angErrNeigh = angErrNeigh+randn()*sigma_measure*pi/LinkNumber;                      % add measurement noise i.i.d. on each agent
        angErrNeigh = angErrNeigh+compassBias(i);                                           % add compass bias (misscalibration) constant for each agent
        
        % compute radial action (u_i,r)
        indices=find(distances > 0 & distances <= MaxSensingRadius);
        v(i,:)= v(i,:) + GainRadial(i) * sum((x(i,:)-xRand(indices,:))./distances(indices) .* RadialInteractionForce(distances(indices), IntFunctionStruct)/InteractionFactor);

        % compute normal action (u_i,n)
        if size(NeigIndices,1) > 0
            v(i,:)= v(i,:) + GainNormal(i) * sum((x(i,:)-xRand(NeigIndices,:))./distances(NeigIndices) * R .* NormalInteractionForce(angErrNeigh, LinkNumber)/InteractionFactor);
        end
        
    end
end