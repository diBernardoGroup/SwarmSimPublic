function [xNeighbours, indices] = getNeighbours(xAgent, x, Rmin, Rmax)
%
%getNeighbours gives the positions and the indeces of the adjacent agents of an agent
%   
%   [xNeighbours, indices] = getNeighbours(xAgent, x, Rmin, Rmax)
%
%   Inputs:
%       xAgent is the position of the agent (1x2 vector)
%       x are the positions of all the agents in the swarm (Nx2 matrix)
%       Rmin and Rmax are the distances defining the adjacency set (scalar)
%
%   Outputs:
%       xNeighbours are the positions of the adjacent agents (matrix)
%       indices are the inideces of the adjacent agents (vector)
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%


    distances=vecnorm(xAgent-x,2,2);

    indices=find((distances<=Rmax) .* (distances>Rmin));
    
    xNeighbours=x(indices,:);
end