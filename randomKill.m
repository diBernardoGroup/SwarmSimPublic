function [xout, sout] = randomKill(x, s, killFactor)
%
%randomKill deletes a fraction of agents selected at random.
%   Used to test the robustness to agents removal.   
%
%   [xout, sout] = randomKill(x, s, killFactor)
%
%   Inputs:
%       x are the positions of the agents (Nx2 matrix)
%       s are the spin of the agents (vector)
%       killFactor is the fraction of agents to delete (scalar [0,1])
%
%   Outputs:
%       xout are the positions of the remaining agents (matrix)
%       sout are the spin of the remaining agents (vector)
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%

    index = randsample(size(x,1),ceil(size(x,1)*(1-killFactor)));
    xout = x(index,:);
    sout = s(index,:);
end

