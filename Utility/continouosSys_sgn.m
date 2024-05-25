function [dx,y] = continouosSys_sgn(t, x, u, theta, alpha, mu, varargin)
    y  = x;
    dx = theta * (mu - x) + sign(x)*u*alpha;
    
%     if x==0
%         dx = 0;
%     else
%         dx = theta * (mu - x) + sign(x)*u*alpha;
%     end

end

