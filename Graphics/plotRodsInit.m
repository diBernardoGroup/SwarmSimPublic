function [p] = plotRodsInit(x,xp,time,windowSize,tickStep,showGrid,thenDelete)
%
%plotSwarm set the correct axis and draws the agents and the links of the swarm.
%
%   [p,p_lines] = plotSwarmInit(x,time,RMin,RMax,Lim,tickStep,showGrid,gradColor,thenDelete)
%
%   Inputs:
%       x           Positions of all the agents                         (NxD matrix)
%       time        Current time instant                                (scalar)       
%       Lim         Size of the window                                  (double = 10)
%       tickStep    Ticks step size                                     (double = Lim/2)
%       showGrid    Display grid                                        (logic = false)
%       gradColor   Use gradient color along the Z axis (3D only)       (logic = false)
%       thenDelete  Delete graphics, used during simulation             (logic = false)
%
%   Outputs:
%       p           Plots of the agents
%
%   See also: plotSwarm, plotTrajectory
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%

arguments
    x           double
    xp          double
    time        double
    windowSize  double {mustBePositive}     = 10
    tickStep    double                      = windowSize/2
    showGrid    logical                     = false
    thenDelete  logical                     = false
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

    [p] = plotRod(x,xp,time,thenDelete);


end

