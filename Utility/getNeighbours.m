function [xNeighbours, indices] = getNeighbours(xAgent, x, Rmin, Rmax)
%
%getNeighbours gives the positions and the indeces of the adjacent agents of an agent
%   
%   [xNeighbours, indices] = getNeighbours(xAgent, x, Rmin, Rmax)
%
%   Inputs:
%       xAgent  Position of the agent                       (1xD vector)
%       x       Positions of all the agents in the swarm    (NxD matrix)
%       Rmin and Rmax are the distances defining the adjacency set (scalar)
%
%   Outputs:
%       xNeighbours Positions of the adjacent agents (matrix)
%       indices     Inideces of the adjacent agents (vector)
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%

arguments
    xAgent      double
    x           double
    Rmin        double {mustBeNonnegative}
    Rmax        double {mustBePositive}
end

    distances=vecnorm(xAgent-x,2,2);

    indices=find((distances<=Rmax) .* (distances>Rmin));
    
    xNeighbours=x(indices,:);
end