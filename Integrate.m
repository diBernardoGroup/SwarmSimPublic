function [xnew, vnew] = Integrate(x, v, Dynamics, deltaT, vMax)
%
%
    
    switch Dynamics.model
    case 'FirstOrder'
        vnew = v;
        noise = Dynamics.sigma * sqrt(deltaT) * randn(size(x));
        
    case 'LevyWalk'
        % select tumbling agents
        tumbling = rand(size(x,1),1)<Dynamics.alpha;
        
        % tumbling agents get a new direction but keep the speed
        velocities=vecnorm(v,2,2);
        new_directions=randn(size(x));
        new_directions=new_directions./vecnorm(new_directions,2,2);
        vnew(tumbling,:)=new_directions(tumbling,:).*velocities(tumbling);
        
        % running agents keep the same velocity
        vnew(~tumbling,:)=v(~tumbling,:);

        noise=zeros(size(x));
        
    otherwise
        error("Dynamics.model is not valid.")
    end
    
    % velocity saturation
    velocities=vecnorm(vnew,2,2);
    indices=find(velocities>vMax);
    vnew(indices,:)= vnew(indices,:) * vMax ./ velocities(indices);
    
    % integration
    xnew = x + vnew.*deltaT + noise;
end

