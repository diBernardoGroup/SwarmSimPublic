function [dx,y] = continouosSys(t, x, u, theta, alpha, mu, varargin)
    y  = x;
    
    if x==0
        dx = 0;
    else
        dx = theta * (mu - x) + u*alpha;
    end
end

