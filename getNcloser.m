function [xSorted] = getNcloser(currPos, x, N)
%
%getNcloser gives the positions of the N closest agents.
%
%   [xSorted] = GetNcloser(currPos, x, N)
%
%   Inputs:
%       currPos is the position to measure from (1x2 vector)
%       x are the positions of all the agents in the swarm (Nx2 matrix)
%       N is the number of agents to return (integer)
%
%   Outputs:
%       xSorted are the positions of the N closest agents, sorted by
%           increasing distance (matrix) 
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%

    xSorted = zeros(N, size(x, 2));
    dist = zeros(size(x,1),1);
    
    for i=1:size(x,1)
        dist(i) = norm(currPos - x(i,:));
        if dist(i) == 0
            dist(i) = inf;
        end
    end
    
    for i=1:N
        [~,ind] = min(dist);
        dist(ind) = inf;
        xSorted(i,:) = x(ind,:);
    end
end

