function [x_inWindow, indices] = getInWindow(x,windowSize)
    inWindow = true(size(x));
    
    for dim = 1:2
        %inWindow(:,dim) = abs(x(:,dim)) <=  windowSize(dim)/2;
        inWindow(:,dim) = x(:,dim) >=  windowSize((dim-1)*2+1) & x(:,dim) <=  windowSize((dim-1)*2+2);
    end
    
    indices = all(inWindow,2);
    x_inWindow = x(indices,:);
end

