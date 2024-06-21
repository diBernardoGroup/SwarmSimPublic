function [p,p_lines] = plotSwarmInit(x,time,RMin,RMax,windowSize,tickStep,showGrid,gradColor,thenDelete, shape, radius, xPrevious)
%
%plotSwarm set the correct axis and draws the agents and the links of the swarm.
%
%   [p,p_lines] = plotSwarmInit(x,time,RMin,RMax,Lim,tickStep,showGrid,gradColor,thenDelete)
%
%   Inputs:
%       x           Positions of all the agents                         (NxD matrix)
%       time        Current time instant                                (scalar)       
%       RMin        Min distance to plot link                           (double)
%       RMax        Min distance to plot link                           (double)
%       Lim         Size of the window                                  (double = 10)
%       tickStep    Ticks step size                                     (double = Lim/2)
%       showGrid    Display grid                                        (logic = false)
%       gradColor   Use gradient color along the Z axis (3D only)       (logic = false)
%       thenDelete  Delete graphics, used during simulation             (logic = false)
%
%   Outputs:
%       p           Plots of the agents
%       p_lines     Plots of the links
%
%   See also: plotSwarm, plotTrajectory
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%

arguments
    x           double
    time        double
    RMin        double {mustBeNonnegative}
    RMax        double {mustBeNonnegative}
    windowSize  double {mustBePositive}     = 10
    tickStep    double                      = windowSize/2
    showGrid    logical                     = false
    gradColor   logical                     = false
    thenDelete  logical                     = false
    shape       string                      = "."
    radius      double {mustBePositive}     = 20
    xPrevious   double                      = x
end
    %figure
    
    if length(windowSize)==1
        windowSize = [windowSize, windowSize];
    end
    
    if length(tickStep)==1
        tickStep = [tickStep, tickStep];
    end
    
    axis('equal',[-windowSize(1)/2 windowSize(1)/2 -windowSize(2)/2 windowSize(2)/2])
    xticks([-windowSize(1)/2:tickStep(1):windowSize(1)/2])
    yticks([-windowSize(2)/2:tickStep(2):windowSize(2)/2])
    xticklabels('')
    yticklabels('')
    zticklabels('')
    box on
    set(gca,'FontSize',14)
    set(gcf,'Position',[100 100 500 500])
    hold on
    if showGrid; grid on; end
    
    if all(size(x) == [2,1])
       x=x'; 
    end

    [p,p_lines] = plotSwarm(x,time, RMin,RMax,thenDelete, ones(size(x,1), 1), gradColor, shape, radius, xPrevious);

end

