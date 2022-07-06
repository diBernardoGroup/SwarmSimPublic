function [p] = drawLines(x,Rmin,Rmax)
%
%drawLines draws the links of the swarm.
%   The figure should be already open and set with the correct axis limits.
%
%   [p] = drawLines(x,Rmin,Rmax)
%
%   Inputs:
%       x are the positions of all the agents (Nx2 matrix)
%       Rmax and Rmin are the distances that define the adjacency set (scalar)
%
%   Outputs:
%       p are the plots of the links
%
%   See also: plotSwarm
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%

 for i=1:length(x)
    for j=i:length(x)
        if norm(x(i,:)-x(j,:))>Rmin && norm(x(i,:)-x(j,:))<Rmax
            p(i,j)=plot([x(i,1) x(j,1)],[x(i,2) x(j,2)], 'k', 'Linewidth', 1);
        end
    end
 end
 
 
 
 end

