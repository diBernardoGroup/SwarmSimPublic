function [p,p_lines,pL] = plotSwarmInit(x,time,RMin,RMax)
%
%plotSwarm set the correct axis and draws the agents and the links of the swarm.
%
%   [p,p_lines,pL] = plotSwarmInit(x,time,RMin,RMax)
%
%   Inputs:
%       x are the positions of all the agents (Nx2 matrix)
%       time is the current time instant (scalar)
%       RMax and RMin are the distances that define the adjacency set (scalar)
%
%   Outputs:
%       p plots of the agents
%       p_lines plots of the links
%       pL plots of the leaders
%
%   See also: plotSwarm
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%

    Max = 10;   % amplitude of the simulation plane
    Min = -Max;

    axis('equal',[Min Max Min Max])
    yticks([-10 -5 0 5 10])
    xticks([-10 -5 0 5 10])
    set(gca,'FontSize',14)
    set(gcf,'Position',[100 100 500 500])
    hold on
    
    if all(size(x) == [2,1])
       x=x'; 
    end

    plotSwarm(x,[],time, RMin,RMax,false, ones(size(x,1), 1));

end

