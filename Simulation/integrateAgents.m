function [xnew, vnew, Dynamics] = integrateAgents(x, v, input, Dynamics, deltaT, envInput)
%
%integrateAgents Integrates agents' dynamics with forward Euler–Maruyama method.
%   You can modify this function to implement your own agents' dynamics.
%   This function is called by Simulator.
%
%   Inputs:
%       x                   Previous positions of the agents    (NxD matrix)
%       v                   Previous velocities of the agents   (NxD matrix)
%       input               Control inputs of the agents        (NxD matrix)
%       Dynamics            Dynamics of the agents              (struct)
%       deltaT              Itegration step                     (positive scalar)
%
%   Outputs:
%       xnew                New positions of the agents         (NxD matrix)
%       vnew                New velocities of the agents        (NxD matrix)
%       Dynamics            Dynamics of the agents              (struct)
%
%   See also: Simulator, Launcher
%
%   Authors:    Andrea Giusti
%   Date:       2023
%

arguments
    x                   double
    v                   double
    input               double
    Dynamics            struct  
    deltaT              double {mustBePositive}
    envInput            double = zeros(size(x,1),1)
end
    
    switch Dynamics.model
    case 'FirstOrder'
        vnew = input;
        noise = Dynamics.sigma * sqrt(deltaT) * randn(size(x));
        
    case 'SecondOrder'
        vnew = v + input * deltaT;
        noise = Dynamics.sigma * sqrt(deltaT) * randn(size(x));
    
    case 'PTW'
        speeds = vecnorm(v,2,2);
        theta = atan2(v(:,2), v(:,1));
        speedsnew = speeds + Dynamics.rateSpeed .* (Dynamics.avgSpeed - speeds) * deltaT + Dynamics.sigmaSpeed * sqrt(deltaT) .* randn(size(x,1),1);
        speedsnew = max(speedsnew, 10e-6);
        Dynamics.omega = Dynamics.omega - Dynamics.rateOmega .* Dynamics.omega * deltaT + Dynamics.sigmaOmega  * sqrt(deltaT) .* randn(size(x,1),1);
        thetanew = mod(theta + pi + Dynamics.omega * deltaT, 2*pi) - pi ;
        vnew = speedsnew .* [cos(thetanew), sin(thetanew)];
        noise = zeros(size(x));
        
    case 'PTWwithInput'
        speeds = vecnorm(v,2,2);
        theta = atan2(v(:,2), v(:,1));
        envInputDotP = max((envInput - Dynamics.oldInput)/deltaT, 0);
        speedsnew = speeds + (Dynamics.rateSpeed .* (Dynamics.avgSpeed-speeds) + Dynamics.gainSpeed .* envInput + Dynamics.gainDerSpeed .* envInputDotP)* deltaT + Dynamics.sigmaSpeed * sqrt(deltaT) .* randn(size(x,1),1);
        speedsnew = max(speedsnew, 10e-6);
        Dynamics.omega = Dynamics.omega + (Dynamics.rateOmega .* (Dynamics.avgOmega-Dynamics.omega) + Dynamics.gainOmega .*envInput + Dynamics.gainDerOmega .* envInputDotP)* deltaT + Dynamics.sigmaOmega  * sqrt(deltaT) .* randn(size(x,1),1);
        thetanew = mod(theta + pi + Dynamics.omega * deltaT, 2*pi) - pi ;
        vnew = speedsnew .* [cos(thetanew), sin(thetanew)];
        noise = zeros(size(x));
        Dynamics.oldInput = envInput;
        
    case 'PTWwithSignedInput'
        speeds = vecnorm(v,2,2);
        theta = atan2(v(:,2), v(:,1));
        b=1/10;
%         envInputDotP = max((power(envInput,b) - power(Dynamics.oldInput,b))/deltaT, 0);
%         envInput = power(envInput,b);
        envInputDotP = max((envInput - Dynamics.oldInput)/deltaT, 0);
        speedsnew = speeds + (Dynamics.rateSpeed .* (Dynamics.avgSpeed - speeds) + Dynamics.gainSpeed .* envInput + Dynamics.gainDerSpeed .* envInputDotP)* deltaT + Dynamics.sigmaSpeed * sqrt(deltaT) .* randn(size(x,1),1);
        speedsnew = max(speedsnew, 10e-6);
        Dynamics.omega = Dynamics.omega + (Dynamics.rateOmega .* (Dynamics.avgOmega-Dynamics.omega) + sign(Dynamics.omega).*(Dynamics.gainOmega .*envInput + Dynamics.gainDerOmega .* envInputDotP))* deltaT + Dynamics.sigmaOmega  * sqrt(deltaT) .* randn(size(x,1),1);
        thetanew = mod(theta + pi + Dynamics.omega * deltaT, 2*pi) - pi ;
        vnew = speedsnew .* [cos(thetanew), sin(thetanew)];
        noise = zeros(size(x));
        Dynamics.oldInput = envInput;
        
    case 'PTWcoupled'
        speeds = vecnorm(v,2,2);
        theta = atan2(v(:,2), v(:,1));
        speedsnew = speeds + Dynamics.rateSpeed * (Dynamics.avgSpeed - speeds) * deltaT + Dynamics.sigmaSpeed * sqrt(deltaT) * randn(size(x,1),1);
        speedsnew = max(speedsnew, 10e-6);
        Dynamics.omega = Dynamics.omega - Dynamics.rateOmega * Dynamics.omega * deltaT + Dynamics.sigmaOmega(speedsnew) * sqrt(deltaT) .* randn(size(x,1),1);
        thetanew = mod(theta + pi + Dynamics.omega * deltaT, 2*pi) - pi ;
        vnew = speedsnew .* [cos(thetanew), sin(thetanew)];
        noise = zeros(size(x));
        
    case 'LevyWalk'
        % select tumbling agents
        tumbling = rand(size(x,1),1)<Dynamics.alpha;
        
        % tumbling agents get a new direction but keep the speed
        if sum(tumbling)
            speeds=vecnorm(v,2,2);
            new_directions=randn(size(x));
            new_directions=new_directions./vecnorm(new_directions,2,2);
            vnew(tumbling,:)=new_directions(tumbling,:).*speeds(tumbling);
        end
        
        % running agents keep the same velocity
        vnew(~tumbling,:)=v(~tumbling,:);
        
        % add input velocity
        vnew = vnew + input;
        
        noise = Dynamics.sigma * sqrt(deltaT) * randn(size(x));
        
    otherwise
        error("Dynamics.model is not valid.")
    end
    
    % velocity saturation
    if isfield(Dynamics,'vMax')
        velocities=vecnorm(vnew,2,2);
        indices=find(velocities>Dynamics.vMax);
        vnew(indices,:)= vnew(indices,:) *Dynamics.vMax ./ velocities(indices);
    end
    
    % integration
    xnew = x + vnew*deltaT + noise;
end

