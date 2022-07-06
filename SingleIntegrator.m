function [xnew] = SingleIntegrator(x, v, deltaT, vMax, sigma)
%
%SingleIntegrator implements the first order dynamics of the agents.
%   Forward Euler with fixed time step integration is used.
%   Saturation of the velocities and random vibrations can be applied.
%
%   [xnew] = SingleIntegrator(x, v, deltaT, vMax, sigma)
%
%   Inputs:
%       x are the positions of the agents (Nx2 matrix)
%       v are the velocities of the agents (Nx2 matrix)
%       deltaT is the integration time step (scalar)
%       vMax is the maximum speed of the agents (scalar)
%       sigma is the standard deviation of noise (scalar)
%
%   Outputs:
%       xnew are the updated positions of the agents (Nx2 matrix)
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%

    % velocity saturation
    velocities=vecnorm(v,2,2);
    indices=find(velocities>vMax);
    v(indices,:)= v(indices,:) * vMax ./ velocities(indices);
    
    % integration
    xnew = x + v.*deltaT + sigma * sqrt(deltaT) * randn(size(x));
end

