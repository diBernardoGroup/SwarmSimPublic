function [p] = plotRods(x,xp,time,thenDelete)
%
%plotSwarm draws the agents and the links of the swarm.
%   The figure should be already open and set with the correct axis using plotSwarmInit.
%
%   [p,p_lines,pL] = plotSwarm(x,xL,time,RMin,RMax,thenDelete, spin, gradColor)
%
%   Inputs:
%       x           Positions of all the agents                         (NxD matrix)   
%       xp          Previous position of all the agents                 (NxD matrix)
%       time        Current time instant                                (scalar)       
%       thenDelete  Delete graphics, used during simulation             (logic = false)
%       spin        Spin of the agents                                  (Nx1 matrix = ones(N,1))
%       gradColor   Use gradient color along the Z axis (3D only)       (logic = false)
%
%   Outputs:
%       p           Plots of the agents
%       p_lines     Plots of the links
%
%   See also: plotSwarmInit, plotTrajectory
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%

arguments
    x           double
    xp          double
    time        double
    thenDelete  logical                     = false
end

warning("plotRods is deprecated, use plotSwarm instead")

title(['t=',num2str(time,'%.2f'),' s'])

if size(x,2) > 3
    x=x';
end


p =gca;
hold on;
for i=1:min([size(xp,1),size(x,1)])
    ang = atan2(x(i,2)-xp(i,2),x(i,1)-xp(i,1));
    plot_singleRod(x(i,:),ang,28,12);
end


if thenDelete
    drawnow
    delete(p)
end

end

