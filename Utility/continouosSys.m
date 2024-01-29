function [dx,y] = continouosSys(t, x, u, theta, alpha, mu, varargin)
    y  = x;
    dx = theta * (mu - x) + u*alpha;
end

