function [v, links, angErr, GainRadial, GainNormal, GainRand] = VFcontroller(x, xL, GainEnv, oldGainRadial, oldGainNormal, oldGainRand, GainLead, GainLeadNormal, RMax, RMin, SensingNumber, InteractionFactor, LinkNumber, deltaT, DeadZoneThresh, IntFunctionStruct, spin, MaxSensingRadius, alpha, beta)
%
%VFcontroller implements the virtual forces controller and computes the
%   velocities (control input) of the agents.
%   This function is called by Simulator.
%
%   [v, links, angErr, GainRadial, GainNormal, GainRand] = ...
%   ... VFcontroller(x, xL, GainEnv, oldGainRadial, oldGainNormal, oldGainRand, GainLead, GainLeadNormal,...
%   ... RMax, RMin, SensingNumber, InteractionFactor, LinkNumber, deltaT, DeadZoneThresh,...
%   ... IntFunctionStruct, spin, MaxSensingRadius, alpha, beta)
%
%   Inputs:
%       x are the positions of the agents (Nx2 matrix)
%       xL are the positions of the leaders (NLx2 matrix)
%       LinkNumber is the desired number of links per agent (scalar)
%       oldGainRadial is the last value of G_radial (vector)
%       oldGainNormal is the last value of G_normal (vector)
%       oldGainRand is the last value of G_random (vector)
%       GainEnv is the value of evironmental gain (scalar)
%       GainLead is the gain of the radial interaction with the leaders (scalar)
%       GainLeadNormal is the gain of the normal interaction with the leaders (scalar)
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
%       GainRand updated values of the random gain (vector)
%
%   See also: Simulator
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%



%% Instantiate Variables
    N=size(x,1);
    v= zeros(N,2);
    links = zeros(N,1);
    angErr = zeros(N,1);
    GainRadial = zeros(N,1);
    GainNormal = zeros(N,1);
    GainRand = zeros(N,1);
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
        
        %% evaluate the phase (only for the hexagonal lattce)
        phase=0;
        if LinkNumber==3
            phase=getPhase(x(i,:), xNeighbours, LinkNumber);
        end
        
        %% links and angular error
        links(i)=size(xNeighbours,1);
        angErr(i)=getAngleError(x(i,:), phase, xNeighbours, LinkNumber);
        
        %% compute gains with adaptation law
        [GainRadial(i), GainNormal(i)]=adaptationLaw(alpha, beta, oldGainRadial(i), oldGainNormal(i), LinkNumber-links(i), angErr(i), deltaT, DeadZoneThresh, LinkNumber);
        
        
        %% compute interactions with neighbours
        % compute the distances from the neighbours (||r_ij||)
        distances=vecnorm(x(i,:)-xRand, 2, 2); 
        
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
        angErrNeigh = getAngularErrNeigh(x(i,:), phase, xNeighbours, LinkNumber);
        
        % compute radial action (u_i,r)
        indices=find(distances > 0 & distances <= MaxSensingRadius);
        v(i,:)= v(i,:) + GainRadial(i) * sum((x(i,:)-xRand(indices,:))./distances(indices) .* RadialInteractionForce(distances(indices), IntFunctionStruct)/InteractionFactor);

        % compute normal action (u_i,n)
        if size(NeigIndices,1) > 0
            v(i,:)= v(i,:) + GainNormal(i) * sum((x(i,:)-xRand(NeigIndices,:))./distances(NeigIndices) * R .* NormalInteractionForce(angErrNeigh, LinkNumber)/InteractionFactor);
        end
        
    end
end