function [p,p_lines] = plotSwarm(x,time,RMin,RMax,thenDelete, spin, gradColor, shape, radius, xPrevious)
%
%plotSwarm draws the agents and the links of the swarm.
%   The figure should be already open and set with the correct axis using plotSwarmInit.
%
%   [p,p_lines,pL] = plotSwarm(x,xL,time,RMin,RMax,thenDelete, spin, gradColor)
%
%   Inputs:
%       x           Positions of all the agents                         (NxD matrix)
%       xL          Positions of the leader agents                      (NLxD matrix)
%       time        Current time instant                                (scalar)       
%       RMin        Min distance to plot link                           (double)
%       RMax        Min distance to plot link                           (double)
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
    time        double                      = 0
    RMin        double {mustBeNonnegative}  = 0
    RMax        double {mustBeNonnegative}  = 0
    thenDelete  logical                     = false
    spin        double                      = ones(size(x,1), 1)
    gradColor   logical                     = true
    shape       string                      = "."
    radius      double {mustBePositive}     = 20
    xPrevious   double                      = x
end

assert(all(size(x)==size(xPrevious)),'x and xPrevious must have the same size!')

if isempty(spin)
        spin = ones(size(x,1), 1);
end

title(['t=',num2str(time,'%.2f'),' s'])
p_lines=drawLines(x,RMin,RMax,gradColor);

spin1 = find(spin==1);
spin0 = find(spin==0);

if size(x,2) > 3
    x=x';
end

if size(x,2) == 2           % 2D plot
    if strcmp(shape,"rod")  % plot rods
        if length(radius)==1; radius=radius*[1,0.5]; end
        ang = atan2(x(:,2)-xPrevious(:,2),x(:,1)-xPrevious(:,1));
        for i=1:size(x,1)
            plot_singleRod(x(i,:),ang(i),radius(1),radius(2));
        end
        ax = gca();
        p1 = ax.Children;
        p2 = [];
    else                    % plot circles
        p1 = plot(x(spin1,1), x(spin1,2),'b.','Marker', shape,'MarkerSize', radius);
        p2 = plot(x(spin0,1), x(spin0,2),'r.','Marker', shape,'MarkerSize', radius);
    end
else                        % 3D plot
    if gradColor
        p1 = scatter3(x(spin1,1), x(spin1,2), x(spin1,3), radius*3, x(spin1,3), 'filled');
        p2 = scatter3(x(spin0,1), x(spin0,2), x(spin0,3), radius*3,'r', 'filled');
        cmap=[0 0 1].*[0.6:0.01:1]';
        colormap(gca,cmap)
    else
        p1 = scatter3(x(spin1,1), x(spin1,2), x(spin1,3),radius*3,'b', 'filled');
        p2 = scatter3(x(spin0,1), x(spin0,2), x(spin0,3),radius*3,'r', 'filled');
    end  
end

p = [p1 ; p2];

if thenDelete
    drawnow
    delete(p)
    delete(p_lines)
end
end

