function [p,p_lines,pL] = plotSwarmInit(x,time,RMin,RMax,Lim,tickStep,showGrid,gradColor)
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

arguments
    x
    time
    RMin
    RMax
    Lim=10
    tickStep=Lim/2
    showGrid=false
    gradColor= false
end

    Max = Lim;   % amplitude of the simulation plane
    Min = -Lim;

    axis('equal',[Min Max Min Max])
    yticks([Min:tickStep:Max])
    xticks([Min:tickStep:Max])
    xticklabels('')
    yticklabels('')
    zticklabels('')
    set(gca,'FontSize',14)
    set(gcf,'Position',[100 100 500 500])
    hold on
    if showGrid; grid on; end
    
    if all(size(x) == [2,1])
       x=x'; 
    end

    plotSwarm(x,[],time, RMin,RMax,false, ones(size(x,1), 1), gradColor);

end

