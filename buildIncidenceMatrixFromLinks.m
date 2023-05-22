function B = buildIncidenceMatrixFromLinks(links,n)
%
%buildIncidenceMatrixFromLinks generates the incidence matrix of a graph starting from the list of links.
%   The list of links can be generated with buildLinks.
%
%   B = buildIncidenceMatrixFromLinks(links,n)
%
%   Inputs:
%       links       List of links as starting and ending vertex     (mx2 matrix)
%       n           Number of vertices                              (integer)
%
%   Outputs:
%       B           Incidence matrix                                (Nxm matrix)
%
%   See also: buildIncidenceMatrix, buildLinks
%
%   Authors:    Andrea Giusti
%   Date:       2023
%

    m=size(links,1);
    B=zeros(n,m);
    
    for k=1:m
        B(links(k,1),k)=1;
        B(links(k,2),k)=-1;
    end
    
    assert(all(size(B)==[n,m]))
end

