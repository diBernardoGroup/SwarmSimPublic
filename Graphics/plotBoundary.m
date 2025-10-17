function [p] = plotBoundary(boundary, color, lineWidth)
    arguments
        boundary  double
        color     = 'k'
        lineWidth = 2
    end
    
    p=plot([-1, 1, 1, -1, -1]*boundary(1)/2, [1, 1, -1, -1, 1]*boundary(2)/2, color=color, lineWidth=lineWidth);
end

