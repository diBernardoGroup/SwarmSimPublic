function [v, links] = VFcontroller(x, GlobalIntFunction, LocalIntFunction, deltaT, InteractionFactor)
%
%VFcontroller implements the virtual forces controller and computes the
%   virtual forces (control input) acting on the agents.
%   This function is called by Simulator.
%
%   [v, links] = VFcontroller(x, GlobalIntFunction, LocalIntFunction, deltaT, InteractionFactor)
%
%   Inputs:
%       x                   Positions of the agents                 (NxD matrix)
%       GlobalIntFunction                                           (struct = struct('function','None'))
%       LocalIntFunction                                            (struct = struct('function','None'))
%       deltaT              Integration time step                   (double = nan)
%       InteractionFactor   Fraction of agents to interact with     (double = 1)
%
%   Outputs:
%       v                   Virtual forces of the agents            (NxD matrix)
%       links               Number of links of each agent           (vector)
%
%   See also: Simulator
%
%   Authors:    Andrea Giusti
%   Date:       2023
%

arguments
    x                   double
    GlobalIntFunction   struct  = struct('function','None')
    LocalIntFunction    struct  = struct('function','None')
    deltaT              double  = nan
    InteractionFactor   double  = 1
end

if ~isfield(LocalIntFunction, 'Rotation')
   LocalIntFunction.Rotation = eye(size(x,2)); 
end

%% Instantiate Variables
    N=size(x,1);
    v= zeros(size(x));
    links = zeros(N,1);
    SensingNumber = GlobalIntFunction.SensingNumber;
    
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
        if isfield(LocalIntFunction,'DistanceRange')
            [xNeighbours, NeigIndices] = getNeighbours(x(i,:), xRand, LocalIntFunction.DistanceRange(1), LocalIntFunction.DistanceRange(2));
        else
            xNeighbours=[];
            NeigIndices=[];
        end
        
        %% links and angular error
        links(i)=size(xNeighbours,1);                
        
        %% compute interactions with neighbours
        % compute the distances from the neighbours (||r_ij||)
        distances=vecnorm(x(i,:)-xRand, 2, 2); 
        
        % Renormalize distances if using Spear'spin algorithm (sqaure or hexagonal lattice)
        if strcmp(GlobalIntFunction.function,'Spears') && (LocalIntFunction.LinkNumber==4 || LocalIntFunction.LinkNumber == 3)
            same_spin = find(spin == spin(i));
            if LinkNumber == 4
                distances(same_spin) = distances(same_spin)/sqrt(2);
            else
                distances(same_spin) = distances(same_spin)/sqrt(3);
            end
        end      
        
        % compute global action (u_i,r)
        if ~strcmp(GlobalIntFunction.function,'None')
            indices=find(distances > 0 & distances <= GlobalIntFunction.MaxSensingRadius);
            v(i,:)= v(i,:) + GlobalIntFunction.Gain * sum((x(i,:)-xRand(indices,:))./distances(indices) .* globalInteractionForce(distances(indices), GlobalIntFunction)/InteractionFactor);
        end
        
        % compute local action (u_i,n)
        if ~strcmp(LocalIntFunction.function,'None')
            if size(NeigIndices,1) > 0
                % compute angular errors (theta_ij^err)
                angErrNeigh = getAngularErrNeigh(x(i,:), 0, xNeighbours, LocalIntFunction.LinkNumber);
                v(i,:)= v(i,:) + LocalIntFunction.Gain * sum((x(i,:)-xRand(NeigIndices,:))./distances(NeigIndices) * LocalIntFunction.Rotation .* localInteractionForce(angErrNeigh, LocalIntFunction.LinkNumber)/InteractionFactor);
            end
        end
    end
end