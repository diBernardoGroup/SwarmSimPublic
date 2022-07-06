function error = getAngleError(xAgent, phase, xNeighbours, LinkNumber)
%
%getAngleError get the average angular error of the agent [normalized].
%
%   error = getAngleError(xAgent, phase, xNeighbours, LinkNumber)
%
%   Inputs:
%       xAgent is the position of the agent (1x2 vector)
%       xNeighbours are the positions of the adjacent agents (matrix)
%       phase is the phase of the agent (scalar)
%       LinkNumber is the desired number of links per agent (scalar)
%
%   Outputs:
%       error is the normalized anglular error of the agent (scalar)
%
%   See also: getAngularErrNeigh, getNeighbours
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%


    refAngle=2*pi/LinkNumber;
    error=0;
    
    if(size(xNeighbours,1)>0)        
        
        bearings= atan2(xAgent(2)-xNeighbours(:,2), xAgent(1)-xNeighbours(:,1))-phase; % angolo rispetto ai vicini (dall'orizzontale)
        closestAngles = round(bearings/refAngle)*refAngle;
        errors=abs(bearings-closestAngles);
        
        error=mean(errors);
        error=error*LinkNumber/pi;
    end
end