function [x_selected, indeces_to_keep, connectivity] = enforceConnectivity(indeces_to_keep, x_lattice, max_iterations)
%
%
%   Authors:    Andrea Giusti
%   Date:       2023


% arguments
%     indeces_to_keep 
%     x_lattice
%     max_iterations double {mustBeInteger, mustBePositive} = 10^3
% end

N = length(indeces_to_keep);
numbers=[1:length(x_lattice)];

% enforce connectivity
counter = 0;
x_selected=x_lattice(indeces_to_keep,:);
B = buildIncidenceMatrix(x_selected, 1.2, false);
connectivity = rank(B)==N-1;
while connectivity==false && counter<max_iterations
    % resample agents
    L=B*B';                         % Laplacian
    A=L-diag(diag(L));              % adjacency matrix
    [bins,binsizes] = conncomp(graph(A));
    [~,largest_concomp]=max(binsizes);
    disconnected_agents = find(bins~=largest_concomp);
    %connected_agents = find(bins==largest_concomp);
    indeces_to_keep(disconnected_agents) = randsample(numbers(all(numbers~=indeces_to_keep)), length(disconnected_agents));
    x_selected=x_lattice(indeces_to_keep,:);
    % check connectivity
    B = buildIncidenceMatrix(x_selected, 1.2, false);
    connectivity = rank(B)==N-1;
    counter = counter+1;
end

% checks
if connectivity==false
    warning("Connectivity not achived")
end
assert(size(x_selected,1)==N)
assert(length(indeces_to_keep) == length(unique(indeces_to_keep)))

%display(['Connectivity achived in iterations ', num2str(counter)])
end

