function y = movmean2(x,k)
%
%movmean2 performs the moving average of a matrix x using a k by k window
%   
%   y = movmean2(x,k)
%
%   Inputs:
%       x is the input matrix (matrix)
%       k is the size of the side of the moving window (integer)
%
%   Outputs:
%       y is the filtered matrix (matrix)
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%

    y=movmean(movmean(x,k,1),k,2);
end

