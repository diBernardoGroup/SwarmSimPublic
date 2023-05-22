function [f] = globalInteractionForce(x, IntFunction)
%
%globalInteractionForce computes the values of radial interaction function for all the neighbours of an agent.
%   You can modify this function to implement your control algorithm.
%
%   [f] = globalInteractionForce(x, IntFunction)
%
%   Inputs:
%       x are the distances between the agent and all its neighbours        (vector)
%       IntFunction describes the interaction function and its parameters   (struct)
%
%   Outputs:
%       f are the values of radial interaction function for all the neighbours of an agent (vector)
%
%   See also: localInteractionForce, VFcontroller
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%

arguments
    x               double
    IntFunction     struct
end

%% evaluate radial interaction function
switch IntFunction.function
    case 'Lennard-Jones'
        a=IntFunction.parameters(1);    
        b=a;
        c=IntFunction.parameters(2);    
        f = a./x.^(2*c)-b./x.^c;        % Lennard-Jones
        f= min(f,1);                    % saturation
        
    case 'Morse'
        %     L=1.1;
        %     F=0.2; %F=1/(exp((L-1)/L));
        F=IntFunction.parameters(1);        %F<1
        L=IntFunction.parameters(2);        %L>1
        S=1; %S=(L-1)/(L*log(1/F));
        f = -F*(exp(-x/L/S))+exp(-x/S);     %Morse
        
    case 'Spears'
        G=IntFunction.parameters(1);
        %K=IntFunction.parameters(2);
        Fmax=IntFunction.parameters(2);
        
        f = G./x.^2; 
        indices=find(x> 1);
        f(indices)= -f(indices);
        indices=find(x> 1.5);
        f(indices)=0;
        f= min(f,Fmax);                % saturation
        f= max(f,-Fmax);               % saturation
        %f=K.*f;
        
    case 'Modified-LJ' %from Torquato2009
        f=(44571*exp(-(249*x)/100))/1000 - (2*exp(-40*(x - 1823/1000).^2).*(80*x - 3646/25))/5 - 589./(10*x.^11) + 60./x.^13;
        f=f/3;
        %f= min(f,1);                % saturation
    
    case 'PowerLaw-FiniteCutoff'
        R=IntFunction.parameters(1);
        Ra=IntFunction.parameters(2);
        f = heaviside(R-x).*(1./x-1/R)*R^2*pi/(Ra-R) - heaviside(Ra-x).*heaviside(x-R).*sin((x-R)*pi/(Ra-R));
    
    otherwise
        error("IntFunction.function must be a valid string ['Lennard-Jones', 'Spears', 'Modified-LJ', 'Morse', 'PowerLaw-FiniteCutoff']")
end
end