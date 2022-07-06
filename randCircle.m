function x = randCircle(N, radius)
%
%randCircle generates N points randomly drawn form a uniform distribution on a circle.
%   Can be used to generate initial positions of the agents.
%
%   x = randCircle(N, radius)
%
%   Inputs:
%       N is the number of points to generate (integer)
%       radius is the radius of the circle (scalar)
%
%   Outputs:
%       x are the random positions (Nx2 matrix)
%
%   See also: perfectLactice
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%

    t = 2*pi*rand(N,1);
    r = radius*sqrt(rand(N,1));
    
    x =  [r.*cos(t), r.*sin(t)];
end

