function [dx,y] = discretizedSys(t, x, u, a, b, c, varargin)
    y  = x;
    dx = a*x + u*b +c;
end

