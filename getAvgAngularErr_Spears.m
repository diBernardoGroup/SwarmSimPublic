function error = getAvgAngularErr_Spears(x, LinkNumber, Rmin, Rmax)
%
%getAvgAngularErr_Spears get the average angular error of the swarm [rad].
%
%   error = getAvgAngularErr_Spears(x, LinkNumber, Rmin, Rmax)
%
%   Inputs:
%       x are the positions of all the agents (Nx2 vector)
%       LinkNumber is the desired number of links per agent (integer)
%       Rmin and Rmax are the distances defining the adjacency set (scalar)
%
%   Outputs:
%       error is the average angular error of the swarm [rad] (scalar)
%
%   See also: getAngleError, getAngularErrNeigh
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%
    
    N=size(x,1);                % number of agents
    refAngle=2*pi/LinkNumber;   % reference angle
    error=0;                    % average angular error of the swarm [rad].
    
    edges=zeros(N*10,1);        % angles of the links of all the agents
    M=0;                        % number of edges, counting both (i,j) and (j,i)
    
    if(N>1)
        for i=1:N
            [xNeighbours] = getNeighbours(x(i,:), x, Rmin, Rmax);
            ph=0;
            if LinkNumber==3;  ph=getPhase(x(i,:), xNeighbours, LinkNumber); end
            bearings = atan2(x(i,2)-xNeighbours(:,2), x(i,1)-xNeighbours(:,1)) - ph; % angles of the links of agent i
            nNeigh=length(bearings);
            edges(M+1:M+nNeigh) = bearings; 
            M=M+nNeigh;
        end
        if M>2
            edges=edges(1:M);
            edgMatrix= edges - edges';
            B = atan2(sin(edgMatrix), cos(edgMatrix));  % matrix of the relative angles between the links
            closestAngles = round(B/refAngle)*refAngle;
            errors=(B-closestAngles);
            
            error = sum(abs(errors),'all')/(M^2-2*M);   % average angular error of the swarm
        end
    end
    
end