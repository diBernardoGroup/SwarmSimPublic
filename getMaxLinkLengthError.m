function maxLinkLenErr = getMaxLinkLengthError(x, R, minSensingRadius, maxSensingRadius)
%
%getMaxLinkLengthError get the normalized maximum absolute error of the length of the links of the swarm (e_d_max). Takes value in [0; max{|Rmax−R|,|Rmin−R|}/R]
%
%   maxLinkLenErr = getMaxLinkLengthError(x, R, minSensingRadius, maxSensingRadius)
%
%   Inputs:
%       x are the positions of all the agent (Nx2 matrix)
%       R is the desired length of the links (scalar)
%       minSensingRadius and maxSensingRadius are the distances defining the adjacency set (scalar)
%
%   Outputs:
%       maxLinkLenErr is the normalized maximum error of the length of the
%           links of the swarm (scalar)
%
%   See also: getAvgLinkLengthError, getAngleError, getNeighbours
%
%   Authors:    Andrea Giusti
%   Date:       2023
%

    D = buildIncidenceMatrix(x, maxSensingRadius);
    r=D'*x;             %realative positions
    d=vecnorm(r,2,2);   %vector of distances
    
    length_errors = abs(d-R);
    
    maxLinkLenErr=max(length_errors)/R;
    
    if isempty(maxLinkLenErr); maxLinkLenErr=nan; end
end