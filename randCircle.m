function x = randCircle(N, radius, D)
%
%randCircle generates N points randomly drawn form a uniform distribution in a disk or a sphere.
%   Can be used to generate initial positions of the agents.
%
%   x = randCircle(N, radius, D)
%
%   Inputs:
%       N is the number of points to generate               (integer)
%       radius is the radius of the circle                  (positive scalar)
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
    radius  double {mustBeNonnegative}
    D       double {mustBeInteger} = 2
end

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

