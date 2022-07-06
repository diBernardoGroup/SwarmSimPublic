function [row, col, costMap] = findOptim(v)
%
%findOptim computes the cost function and finds its minimum.
%   This function is called by BruteForceTuning.
%
%   [row, col, costMap] = findOptim(v)
%
%   Inputs:
%       v are the steady state values of the metrics (3D matrix)
%
%   Outputs:
%       row  and col are the indices of the minimum (integer)
%       costMap stores the avlues of the cost function (2D matrix)
%
%   See also: BruteForceTuning
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%

    % 2 norm
    L=vecnorm(v,2,3);
    [~, indx]=min(L(:));
    [row, col]=ind2sub(size(L),indx);

    
%     % inifinity norm (maximum)
%     L=vecnorm(v,inf,3);
%     %L(L>1)=1;
%     [~, indx]=min(L(:));
%     [row, col]=ind2sub(size(L),indx);
    
    costMap=L;
end

