function [x, v, out_agents] = boundaryInteraction(x, v, boundary)
    out_agents = false(size(x));
    for dim = 1:length(boundary)
        out_agents(:,dim) = abs(x(:,dim)) >  boundary(dim)/2;
    end
    
    for dim = 1:length(boundary)
        x(out_agents(:,dim),dim) = min(x(out_agents(:,dim),dim), boundary(dim)/2);
        x(out_agents(:,dim),dim) = max(x(out_agents(:,dim),dim), -boundary(dim)/2);
        
        %v(out_agents(:,dim),dim) = -x(out_agents(:,dim))/10000; % boundary repulsion
        v(out_agents(:,dim),dim) = 0;
    end
    
    
%     for dim = 1:length(boundary)
%         x(:,dim) = min(x(:,dim), boundary(dim)/2);
%         x(:,dim) = max(x(:,dim), -boundary(dim)/2);
%     end
end

