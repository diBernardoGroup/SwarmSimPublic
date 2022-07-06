function phase = getPhase(xAgent, xNeighbours, LinkNumber)
%
%getPhase gives the angular phase of the agent [rad].
%   Used only for hexagonal lattice.
%
%   phase = getPhase(xAgent, xNeighbours, LinkNumber)
%
%   Inputs:
%       xAgent is the position of the agent (1x2 vector)
%       xNeighbours are the positions of the adjacent agents (matrix)
%       LinkNumber is the desired number of links per agent (scalar)
%
%   Outputs:
%       phase is the phase of the agent (scalar)
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%

    % role evaluation
    phase1=0;
    phase2=pi*LinkNumber;
    cost1=0; cost2=0;
    
    angErrNeig1 = getAngularErrNeigh(xAgent, phase1, xNeighbours, LinkNumber);
    angErrNeig2 = getAngularErrNeigh(xAgent, phase2, xNeighbours, LinkNumber);
    
    cost1 = sum(angErrNeig1.^2);
    cost2 = sum(angErrNeig2.^2);

    % phase assignment
    if cost1<cost2; phase=phase1;
    else; phase=phase2;
    end

end

