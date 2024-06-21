function [x_inWindow, indices] = getInWindow(x,windowSize)
    inWindow = true(size(x));
    
    for dim = 1:length(windowSize)
        inWindow(:,dim) = abs(x(:,dim)) <=  windowSize(dim)/2;
    end
    
    indices = all(inWindow,2);
    x_inWindow = x(indices,:);
end

