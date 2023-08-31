function errors = getAngularErrNeigh(xAgent, phase, xNeighbours, LinkNumber)
%
%getAngularErrNeigh get the angular errors of the links of the agent [rad].
%
%   errors = getAngularErrNeigh(xAgent, phase, xNeighbours, LinkNumber)
%
%   Inputs:
%       xAgent is the position of the agent (1x2 vector)
%       xNeighbours are the positions of the adjacent agents (matrix)
%       phase is the phase of the agent (scalar)
%       LinkNumber is the desired number of links per agent (integer)
%
%   Outputs:
%       errors are the anglular errors of the links (vector)
%
%   See also: getAngleError, getNeighbours
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%


    refAngle=2*pi/LinkNumber;
    errors=zeros(size(xNeighbours,1));
    
    if(size(xNeighbours,1)>0)
        bearings= atan2(xAgent(2)-xNeighbours(:,2), xAgent(1)-xNeighbours(:,1))-phase; % angolo dei vicini (dall'orizzontale)
        closestAngles = round(bearings/refAngle)*refAngle;
        errors=(bearings-closestAngles);
    end
end