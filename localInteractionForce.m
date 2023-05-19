function [f] = localInteractionForce(theta,LinkNumber)
%
%localInteractionForce computes the values of normal interaction function for all the links of an agent.
%
%   [f] = localInteractionForce(theta,LinkNumber)
%
%   Inputs:
%       theta are the angular errors of all the links of the agent (vector)
%       LinkNumber is the desired number of links per agent (scalar)
%
%   Outputs:
%       f are the values of normal interaction function for all the links of an agent (vector)
%
%   See also: globalInteractionForce, VFcontroller
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%


    %f=0;                                       % null interaction
    %f=-sin((theta+phase)*LinkNumber)/d^2;       % sin
    %f=-min(max(tan(theta*LinkNumber/2),-1),1);  % saturated tangent
    f= -theta*LinkNumber/pi;                     % linear from 1 to -1 between -pi/L and pi/L
end