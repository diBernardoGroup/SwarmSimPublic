function B = buildIncidenceMatrix(x, Rmax, oriented)
%
%buildIncidenceMatrix generates the incidence matrix of the proximity graph.
%
%   B = buildIncidenceMatrix(x, Rmax, oriented)
%
%   Inputs:
%       x           Positions of the verices                            (Nx: matrix)
%       Rmax        Max length of a link                                (scalar)
%       oriented    If true each link is taken with both directions     (logical = true)
%
%   Outputs:
%       B           Incidence matrix                                    (Nxm matrix)
%
%   See also: buildIncidenceMatrixFromLinks
%
%   Authors:    Andrea Giusti
%   Date:       2023
%

arguments
    x           double
    Rmax        double {mustBePositive}
    oriented    logical = true
end

n=size(x,1);
m=0;

B=zeros(n,0);

% for each agent
for i=1:n 
    xAgent=x(i,:);
    [~, indices] = getNeighbours(xAgent, x, 10e-6, Rmax);
    
    % for each of its neighbours
    for j=[indices]'
        if oriented || j>i
            m=m+1;
            B(i,m)=1;
            B(j,m)=-1;
        end
    end
    
end

assert(size(B,1)==n, "D is "+size(B,1)+"x"+size(B,2)+" but len(X)="+size(x,1))
assert(all(sum(B)==zeros(1,m)),"Columns of D do not sum to zero.")

end

