function [f] = localInteractionForce(xAgent, xNeighbours, LocalIntFunction)
%
%localInteractionForce computes the magnitude of the local interaction forces generated by the neighbours of an agent.
%   You can modify this function to implement your control algorithm.
%
%   [f] = localInteractionForce(xAgent, xNeighbours, LocalIntFunction)
%
%   Inputs:
%       xAgent is the position of the agent                             (vector)
%       xNeighbours are the positions of the neighbours of the agent    (matrix)
%       LocalIntFunction describes shorts distance interaction          (struct)
%
%   Outputs:
%       f are the magnitude of the local interaction forces applied by the neighbours of an agent (vector)
%
%   See also: globalInteractionForce, VFcontroller
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%

arguments
    xAgent      (1,:)   double
    xNeighbours         double
    LocalIntFunction    struct
end

L=LocalIntFunction.LinkNumber;
theta = getAngularErrNeigh(xAgent, 0, xNeighbours, L);

switch LocalIntFunction.function
    case 'Sin'
        f=-sin((theta)*L);                  % sin
        
    case 'Tan'
        f=-min(max(tan(theta*L/2),-1),1);   % saturated tangent
        
    case 'Linear'
        f= -theta*L/pi;                     % linear from 1 to -1 between -pi/L and pi/L
        
    otherwise
        error("LocalIntFunction.function must be a valid string ['Sin', 'Tan', 'Linear']")
end

if size(f,1)<size(f,2); f=f'; end

end
