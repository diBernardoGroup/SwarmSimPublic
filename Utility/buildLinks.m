function links = buildLinks(x, Rmax, oriented)
%
%buildLinks generates the list of links of the proximity graph.
%
%   links = buildLinks(x, Rmax, oriented)
%
%   Inputs:
%       x           Positions of the verices                            (Nx: matrix)
%       Rmax        Max length of a link                                (scalar)
%       oriented    If true each link is taken with both directions     (logical = true)
%
%   Outputs:
%       links       List of links                                       (mx2 matrix)
%
%   See also: buildIncidenceMatrixFromLinks
%
%   Authors:    Andrea Giusti
%   Date:       2023
%

arguments
    x
    Rmax
    oriented = true
end

    n=size(x,1);
    m=0;
    
    links=zeros(0,2);
    
    for i=1:n
        xAgent=x(i,:);
        [~, indices] = getNeighbours(xAgent, x, 10e-6, Rmax);
       
        for j=[indices]'
            if oriented || j>i
                m=m+1;
                links(m,:)=[i,j];
            end
        end
        
    end 
end

