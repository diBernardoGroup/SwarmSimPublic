function [theta, alpha, mu, sigma] = dt2ct_parameters(a,b,c,s, deltaT)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
  %compute continuous time parameters
            mu = real(c/(1-a));
            alpha = real(inv(a-1)*log(a)/deltaT*b);
            if a ~= 0
                theta = real(-1/deltaT * log(a));
                sigma = real(s * sqrt(-2*log(a) / (1-a^2) / deltaT));
            else
                theta = 0;
                sigma = 0;
            end
end

