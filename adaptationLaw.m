function [newGainRadial, newGainNormal] = adaptationLaw(alpha, beta, GainRadial, GainNormal, linkError, NormalisedAngleError, deltaT, DeadZoneThresh, LinkNumber)
%
%adaptationLaw computes the updated control gains of an agent.
%   Inputs:
%       alpha and beta are the adaptation gains (scalar)
%       GainRadial is the previous value of the radial gain (scalar)
%       GainNormal is the previous value of the normal gain (scalar)
%       linkError is the difference betweeen LinkNumber and the number of
%           links of the agent (scalar)
%       NormalisedAngleError is the normalised angular error of the agent
%           (scalar) (theta_i^err)
%       deltaT integration time step (scalar)
%       DeadZoneThresh amplitute of the adaptation deadhzone (scalar)
%       LinkNumber is the desired number of links per agent (scalar)
%
%   Outputs:
%       newGainRadial is the updated value of the radial gain (scalar)
%       newGainNormal is the updated value of the normal gain (scalar)
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%


    LinkTresh = 0;
    
    newGainRadial = GainRadial;
    newGainNormal = GainNormal;
    
    if NormalisedAngleError > DeadZoneThresh
        %1
        % newGainNormal = GainNormal + alpha * (NormalisedAngleError - DeadZoneThresh) * deltaT;
        
        %2
        newGainNormal = GainNormal + (alpha * (NormalisedAngleError - DeadZoneThresh) - beta * linkError) * deltaT;
        
        %3
        % newGainNormal = GainNormal + (alpha * (NormalisedAngleError - DeadZoneThresh)/(linkError+1)) * deltaT;
        
        %4
        % if linkError <= LinkTresh
        %    newGainNormal = GainNormal + alpha * (NormalisedAngleError - DeadZoneThresh) * deltaT;
        % else
        %    newGainNormal = GainNormal;
        % end
        
        %5
        % if linkError <= LinkTresh
        %    newGainNormal = GainNormal + alpha * (NormalisedAngleError - DeadZoneThresh) * deltaT;
        % else
        %    newGainNormal = GainNormal - beta * linkError * deltaT;
        % end
    end
    
    newGainRadial = max(0, newGainRadial);
    newGainNormal = max(0, newGainNormal);

end

