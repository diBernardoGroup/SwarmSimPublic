function M = buildRigidityMatrix(x, B)
%
%buildRigidityMatrix generates the rigidity matrix of the graph.
%   If rank(M)==Dn-D(D+1)/2 the graph is infinitesimally rigid.
%
%   M = buildRigidityMatrix(x, B)
%   
%   Inputs:
%       x Positions of the vertices (NxD matrix)
%       B Incidence matrix          (Nxm matrix)
%
%   Outputs:
%       M Rigidity matrix           (mx(DN) matrix)
%
%   See also: buildIncidenceMatrix
%
%   Authors:    Andrea Giusti
%   Date:       2022
%

arguments
    x   double
    B   double
end

    n=size(B,1);
    m=size(B,2);
    D=size(x,2);
    
    R=B'*x;
    M=kron(B', ones(1,D)) .* kron(ones(1,n), R);
 
    assert(all(size(M)==[m, D*n]),"size(M)="+size(M));
end

