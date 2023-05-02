function x = randCircle(N, radius, D)
%
%randCircle generates N points randomly drawn form a uniform distribution on a circle.
%   Can be used to generate initial positions of the agents.
%
%   x = randCircle(N, radius, D)
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

azimuth = 2*pi*rand(N,1);
r = radius*nthroot(rand(N,1),D);

if D==2
    x = [r.*cos(azimuth), r.*sin(azimuth)];
else
    rvals = 2*rand(N,1)-1;
    elevation = asin(rvals);
    [x,y,z] = sph2cart(azimuth,elevation,r);
    
    x = [x,y,z];
end
end

