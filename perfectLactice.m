function x_selected = perfectLactice(N, L, D, enf_connectivity, enf_rigidity, lattice_size, max_iterations)
%
%perfectLactice generates N points arranged on a lattice.
%   Can be used to generate initial positions of the agents.
%
%   x = perfectLactice(N, L, D, enf_connectivity, enf_rigidity, lattice_size, max_iterations)
%
%   Inputs:
%       N                   Number of points to generate                                        (integer)
%       L                   Number of links per agent                                           (integer)
%           If D=2 valid values are 3=hexagonal lattice, 4=sqaure lattice, 6=triangular lattice
%           If D=3 valid values are 6=cubic lattice, 12=thetradic-octaedric lattice
%       D                   Dimension of the space. Must be 2 or 3.                             (integer)
%       enf_connectivity    Ensure the resulting swarm graph is connected                       (logical = true)
%       enf_rigidity        Ensure the resulting swarm graph is rigid                           (logical = true)
%       lattice_size        Number of points in the lattice, must be larger or equal than N     (integer = ceil(sqrt(N))^2)
%       max_iterations      Maximum number of iterations                                        (positive integer = 10^3)
%   
%   Outputs:
%       x_selected          Positions of the selected lattice points                            (NxD matrix)
%
%   See also: randCircle
%
%   Authors:    Andrea Giusti
%   Date:       2023
%

arguments
    N double {mustBeInteger, mustBePositive}
    L double {mustBeInteger, mustBePositive, mustBeMember(L,[3,4,6,12])}
    D double {mustBeInteger, mustBePositive, mustBeMember(D,[2,3])} = 2
    enf_connectivity logical = true
    enf_rigidity logical = true
    lattice_size double {mustBeInteger, mustBePositive, mustBeGreaterThanOrEqual(lattice_size,N)} = ceil(nthroot(N,D))^D
    max_iterations double {mustBeInteger, mustBePositive} = 10^3
end

if enf_rigidity && ~enf_connectivity
    enf_connectivity=true;
    warning('A graph cannot be rigid and disconeected. enf_connectivity has been set to true.')
end

if D==2
    assert(ismember(L,[3,4,6]),'If D=2 then L must belong to [3,4,6].')
    assert(floor(sqrt(lattice_size))==ceil(sqrt(lattice_size)),'If D=2 then lattice_size must be a square number.')
    
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
    
elseif D==3
    assert(ismember(L,[6,12]),'If D=3 then L must belong to [6,12].')
    assert(floor(nthroot(lattice_size,3))==ceil(nthroot(lattice_size,3)),'If D=3 then lattice_size must be a cubic number.')
    
    l=nthroot(lattice_size,3);
    l = ceil( l );
    
    % create the lattice
    switch L
        case 6
            
            for i=1:l
                for j=1:l
                    for k=1:l
                        index=sub2ind([l,l,l],i,j,k);
                        x(index,:)=[(i-1), (j-1), (k-1)];
                    end
                end
            end
            
        case 12
            
            for i=1:l
                for j=1:l
                    for k=1:l
                        index=sub2ind([l,l,l],i,j,k);
                        x(index,:)=[(i-1)+mod(j,2)/2+mod(k,2)/2, (j-1)*sqrt(3)/2+mod(k,2)/sqrt(12), (k-1)*sqrt(2/3)];
                    end
                end
            end
            
    end
end

% sample the required number of agents
indeces_to_keep = randsample(length(x), N);
x_selected=x(indeces_to_keep,:);

numbers=[1:length(x)];

% enforce connectivity
if enf_connectivity
    [x_selected, indeces_to_keep, connectivity] = enforceConnectivity(indeces_to_keep, x, max_iterations);
end

% enforce rigidity
counter = 0;
B = buildIncidenceMatrix(x_selected, 1.2, false);
M = buildRigidityMatrix(x_selected, B);
rigidity = rank(M)==D*N-D*(D+1)/2;
while enf_rigidity && rigidity==false && counter<max_iterations
    % resample agents with fewer links
    [~,agents_to_resample]=mink(sum(abs(B')),ceil(N/5));
    agents_to_keep = setdiff(indeces_to_keep, indeces_to_keep(agents_to_resample));
    indeces_to_keep(agents_to_resample) = randsample(setdiff(numbers,agents_to_keep), length(agents_to_resample));
    x_selected=x(indeces_to_keep,:);
    B = buildIncidenceMatrix(x_selected, 1.2, false);
    
    % enforce connectivity and at least D links per agent
    counter2=0;
    while (~connectivity || min(sum(abs(B')))<D) && counter2<ceil(max_iterations/100)+10
        [x_selected, indeces_to_keep, connectivity] = enforceConnectivity(indeces_to_keep, x, ceil(max_iterations/100)+10);
        B = buildIncidenceMatrix(x_selected, 1.2, false);
        peripheric_agents = find(sum(abs(B'))<min([D,sum(abs(B'))+1]));
        central_agents = setdiff(indeces_to_keep, indeces_to_keep(peripheric_agents));
        indeces_to_keep(peripheric_agents) = randsample(setdiff(numbers,central_agents), length(peripheric_agents));
        [x_selected, indeces_to_keep, connectivity] = enforceConnectivity(indeces_to_keep, x, ceil(max_iterations/100)+10);
        B = buildIncidenceMatrix(x_selected, 1.2, false);
        counter2 = counter2+1;
    end
    % check rigidity
    M = buildRigidityMatrix(x_selected, B);
    rigidity = rank(M)==D*N-D*(D+1)/2;
    counter = counter+1;
    % display(['Minimal number of links achived in iterations ', num2str(counter2)])
end

% enforce connectivity
if enf_connectivity && ~connectivity
    [x_selected, indeces_to_keep, connectivity] = enforceConnectivity(indeces_to_keep, x, max_iterations);
end

% recenter the lattice
x_selected=x_selected-mean(x_selected);

% checks
if enf_connectivity && connectivity==false
    warning("Connectivity not achived")
end
if enf_rigidity && rigidity==false
    warning("Rigidity not achived")
end

assert(all(size(x_selected)==[N,D]))
assert(length(indeces_to_keep) == length(unique(indeces_to_keep)))

%display(['Rigidity achived in iterations ', num2str(counter)])
end

