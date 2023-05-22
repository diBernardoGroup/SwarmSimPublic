function avgLinkLenErr = getAvgLinkLengthError(x, R, minSensingRadius, maxSensingRadius)
%
%getAvgLinkLengthError get the normalized average error of the length of the links of the swarm (e_d). Takes value in [0; max{|Rmax−R|,|Rmin−R|}/R]
%
%   avgLinkLenErr = getAvgLinkLengthError(x, R, minSensingRadius, maxSensingRadius)
%
%   Inputs:
%       x are the positions of all the agent (Nx2 matrix)
%       R is the desired length of the links (scalar)
%       minSensingRadius and maxSensingRadius are the distances defining the adjacency set (scalar)
%
%   Outputs:
%       avgLinkLenErr is the normalized average error of the length of the
%           links of the swarm (scalar)
%
%   See also: getAngleError, getNeighbours
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%

    N=size(x,1);
    
    avgErr=zeros(N,1);
    
    if(N>1)
        for i=1:N
            [xNeighbours] = getNeighbours(x(i,:), x, minSensingRadius, maxSensingRadius);
            if(~isempty(xNeighbours))
                distances = vecnorm(x(i,:)-xNeighbours, 2, 2); 
                errors=abs(distances-R);
                avgErr(i)=mean(errors);
            end
        end
    end
    
    avgLinkLenErr=mean(avgErr)/R;
    
end