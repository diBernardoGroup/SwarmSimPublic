function [p,p_lines,pL] = plotSwarm(x,xL,time,RMin,RMax,thenDelete, spin)
%
%plotSwarm draws the agents and the links of the swarm.
%   The figure should be already open and set with the correct axis limits.
%
%   [p,p_lines,pL] = plotSwarm(x,xL,time,RMin,RMax,thenDelete, spin)
%
%   Inputs:
%       x are the positions of all the agents (Nx2 matrix)
%       xL are the positions of the leaders (NLx2 matrix)
%       time is the current time instant (scalar)
%       RMax and RMin are the distances that define the adjacency set (scalar)
%       thenDelete must be true to make animations (bool)
%       spin are the spin of the agents (vector of bools)
%
%   Outputs:
%       p plots of the agents
%       p_lines plots of the links
%       pL plots of the leaders
%
%   See also: drawLines
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%

    title("t="+time+" s")
    p_lines=drawLines(x,RMin,RMax);
    
    spin1 = find(spin==1);
    p1 = plot(x(spin1,1), x(spin1,2),'b.','MarkerSize', 20);
    spin0 = find(spin==0);
    p2 = plot(x(spin0,1), x(spin0,2),'r.','MarkerSize', 20);
    p = [p1 ; p2];
    
    if(size(xL,1)>0); pL = plot(xL(:,1), xL(:,2),'r.','MarkerSize', 20); end

    if(thenDelete)
        drawnow
        delete(p)
        delete(p_lines)
        if(size(xL,1)>0); delete(pL); end
    end
end

