function x = perfectLactice(N, L)
%
%perfectLactice generates N points arranged on a lattice.
%   Can be used to generate initial positions of the agents.
%
%   x = perfectLactice(N, L)
%
%   Inputs:
%       N is the number of points to generate (integer)
%       L is the number of links per agent (integer)
%           valid values are 3=hexagonal lattice, 4=sqaure lattice, 6=triangular lattice
%           if L==4 then N must be a square number.
%
%   Outputs:
%       x are the positions of the lattice points (Nx2 matrix)
%
%   See also: randCircle
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%

    x=nan(N,2);

    switch L
        case 4
            l=sqrt(N);

            if floor(l)==ceil(l)
                for i=1:l
                    for j=1:l
                        index=(i-1)*l+j;
                        x(index,:)=[(i-1), (j-1)];
                    end
                end
            else
                display('!!!!!! N must be a square number. !!!!!!');
            end

        case 6
            l=sqrt(N);

            for i=1:l
                for j=1:l
                    index=(i-1)*l+j;
                    x(index,:)=[(i-1)+mod(j,2)/2, (j-1)*sqrt(3)/2];
                end
            end
            
        case 3
            l=sqrt(N);
            
            for i=1:l
                for j=1:l
                    index=(i-1)*l+j;
                    x(index,:)=[(i-1)/2+floor((i-1)/2)/2, (j-1)*sqrt(3)+mod(floor((i)/2),2)*sqrt(3)/2];
                end
            end
    end

    x=x-mean(x);
end

