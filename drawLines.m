function [p] = drawLines(x, RMin, RMax, gradColor)
%
%drawLines draws the links of the swarm.
%   The figure should be already open and set with the correct axis using plotSwarmInit.
%
%   [p] = drawLines(x, RMin, RMax, gradColor)
%
%   Inputs:
%       x           Positions of all the agents                         (NxD matrix)
%       RMin        Min distance to plot link                           (double)
%       RMax        Min distance to plot link                           (double)
%       gradColor   Use gradient color along the Z axis (3D only)       (logic = false)
%
%   Outputs:
%       p           Plots of the links
%
%   See also: plotSwarm, plotSwarmInit
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%

arguments
    x           double
    RMin        double {mustBeNonnegative}
    RMax        double {mustBeNonnegative}
    gradColor   logical                     = false
end

p=[];

for i=1:size(x,1)
    for j=i:size(x,1)
        if norm(x(i,:)-x(j,:))>RMin && norm(x(i,:)-x(j,:))<RMax
            if size(x,2) == 2
                p(i,j)=plot([x(i,1) x(j,1)],[x(i,2) x(j,2)], 'k', 'Linewidth', 1);
            else
                if gradColor
                    p(i,j)=surface([[x(i,1) x(j,1)];[x(i,1) x(j,1)]],[[x(i,2) x(j,2)];[x(i,2) x(j,2)]],[[x(i,3) x(j,3)];[x(i,3) x(j,3)]],[[x(i,3) x(j,3)];[x(i,3) x(j,3)]],'facecol','no','edgecol','interp','linew',1);
                else
                    p(i,j)=plot3([x(i,1) x(j,1)],[x(i,2) x(j,2)],[x(i,3) x(j,3)], 'k', 'Linewidth', 1);
                end
            end
        end
    end
 end
 
 
 
 end

