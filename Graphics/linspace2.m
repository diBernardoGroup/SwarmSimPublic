function [matrix] = linspace2(start, stop, npoints)
    matrix = nan(length(start), npoints);

    for i=1:length(start)
        matrix(i,:)= linspace(start(i),stop(i),npoints);
    end
end

