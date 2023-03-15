function x_selected = perfectLactice(N, L, enf_connectivity, enf_rigidity, lattice_size, max_iterations)
%
%perfectLactice generates N points arranged on a lattice.
%   Can be used to generate initial positions of the agents.
%
%   x = perfectLactice(N, L, enf_connectivity, enf_rigidity, lattice_size)
%
%   Inputs:
%       N is the number of points to generate. N must be a square number. (integer)
%       L is the number of links per agent (integer)
%           valid values are 3=hexagonal lattice, 4=sqaure lattice, 6=triangular lattice
%       enf_connectivity ensures the resulting swarm graph is connected (logical = true)
%       enf_rigidity ensures the resulting swarm graph is rigid (logical = true)
%       lattice_size number of points in the lattice, must be square and larger or equal than N (integer = ceil(sqrt(N))^2)
%   Outputs:
%       x_selected are the positions of the lattice points (Nx2 matrix)
%
%   See also: randCircle
%
%   Authors:    Andrea Giusti
%   Date:       2023
%
    arguments
        N double {mustBeInteger, mustBePositive}
        L double {mustBeInteger, mustBePositive, mustBeMember(L,[3,4,6])}
        enf_connectivity logical = true
        enf_rigidity logical = true
        lattice_size double {mustBeInteger, mustBePositive, mustBeGreaterThanOrEqual(lattice_size,N)} = ceil(sqrt(N))^2
        max_iterations double {mustBeInteger, mustBePositive} = 10^3
    end
    
    if enf_rigidity; enf_connectivity=true; end
    
    assert(floor(sqrt(lattice_size))==ceil(sqrt(lattice_size)),'lattice_size must be a square number.')

    l=sqrt(lattice_size);
    l = ceil( l );
    
    % create the lattice
    switch L
        case 4

            for i=1:l
                for j=1:l
                    index=(i-1)*l+j;
                    x(index,:)=[(i-1), (j-1)];
                end
            end

        case 6

            for i=1:l
                for j=1:l
                    index=(i-1)*l+j;
                    x(index,:)=[(i-1)+mod(j,2)/2, (j-1)*sqrt(3)/2];
                end
            end
            
        case 3
            
            for i=1:l
                for j=1:l
                    index=(i-1)*l+j;
                    x(index,:)=[(i-1)/2+floor((i-1)/2)/2, (j-1)*sqrt(3)+mod(floor((i)/2),2)*sqrt(3)/2];
                end
            end
    end
    
    % sample the required number of agents
    if length(x) > N
        indeces_to_keep = randsample(length(x), N);
        x_selected=x(indeces_to_keep,:);
    else
        x_selected=x;
    end
    
    % enforce rigidity
    counter = 0;
    D = buildIncidenceMatrix(x_selected, 1.2);
    M = buildRigidityMatrix(x_selected, D);
    rigidity = rank(M)==2*N-3;
    while enf_rigidity && rigidity==false && counter<max_iterations
        % resample agents
        indeces_to_keep = randsample(length(x), N);
        x_selected=x(indeces_to_keep,:);
        % check rigidity
        D = buildIncidenceMatrix(x_selected, 1.2);
        M = buildRigidityMatrix(x_selected, D);
        rigidity = rank(M)==2*N-3;
        counter = counter+1;
    end
    
    % enforce connectivity
    counter = 0;
    D = buildIncidenceMatrix(x_selected, 1.2);
    connectivity = rank(D)==N-1;
    while enf_connectivity && connectivity==false && counter<max_iterations
        % resample agents
        indeces_to_keep = randsample(length(x), N);
        x_selected=x(indeces_to_keep,:);
        % check connectivity
        D = buildIncidenceMatrix(x_selected, 1.2);
        connectivity = rank(D)==N-1;
        counter = counter+1;
    end 
    
    % checks
    if enf_connectivity && connectivity==false
        warning("Connectivity not achived")
    end
    if enf_rigidity && rigidity==false
        warning("Rigidity not achived")
    end
    
    % recenter the lattice
    x_selected=x_selected-mean(x_selected);
    
    assert(size(x_selected,1)==N)
end

