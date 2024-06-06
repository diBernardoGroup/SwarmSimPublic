function x = randRect(N, sizes, D)
%
%randRect generates N points randomly drawn form a uniform distribution onto a rectangle or a cube.
%   Can be used to generate initial positions of the agents.
%
%   x = randRect(N, sizes, D)
%
%   Inputs:
%       N is the number of points to generate               (integer)
%       sizes is the radius of the circle                  (positive scalar)
%       D is the dimension of the space. Must be 2 or 3.    (integer = 2)
%
%   Outputs:
%       x are the random positions                          (NxD matrix)
%
%   See also: perfectLactice
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%

arguments
    N       double {mustBeInteger}
    sizes   double {mustBeNonnegative}
    D       double {mustBeInteger} = 2
end

    x = (rand(N,length(sizes))-0.5) .* sizes;
end

