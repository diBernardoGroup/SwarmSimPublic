function [T_r, T_theta, T_L, success, k1, k2] = settlingTimes(e_theta,e_L, regularity_thresh, compactness_thresh, TSample)
%
%settlingTimes computes the rising times given the time-series of the metrics.
%
%   [T_r, T_theta, T_L, success, k1, k2] = settlingTimes(e_theta,e_L, regularity_thresh, compactness_thresh, TSample)
%
%   Inputs:
%       e_theta is the regularity metrics (vector)
%       e_L is the compactness metrics (vector)
%       regularity_thresh is the threshold value for regularity metrics (e^*_theta) (scalar)
%       compactness_thresh is the threshold value for compactness metrics (e^*_L) (scalar)
%       TSample are the sampling time instants (vector)
%
%   Outputs:
%       T_r is the overall rising time (scalar)
%       T_theta is the rising time of the regularity metrics (scalar)
%       T_L is the rising time of the compactness metrics (scalar)
%       success is true if both metrics get below the respecive threshold (bool)
%       k1 is the index of the time instant corresponding to T_theta (integer)
%       k2 is the index of the time instant corresponding to T_L (integer)
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%

    success=true; 
    
    %T_theta
    k1=find(movmean(e_theta,3) >= regularity_thresh, 1, 'last') + 1;
    if(k1 > length(e_theta))
        success=false; 
        k1=NaN;
        T_theta=NaN;
    else
        T_theta=TSample(k1);
    end

    %T_L
    k2=find(movmean(e_L,3) >= compactness_thresh, 1, 'last') + 1;
    if(k2 > length(e_L))
        success=false; 
        k2=NaN;
        T_L=NaN;
    else
        T_L=TSample(k2);
    end
    
    T_r=max(T_L,T_theta,'includenan');
    
    %T_ALN deprecated
    %     errMeanALN = mean(avgLinkNumber(1,length(avgLinkNumber)-10:end));
    %     k3=find(avgLinkNumber(1,:) < errMeanALN*1.2, 1);
    %     if(isempty(k3)); k3=length(TSample);end
    %     T_ALN=TSample(k3);
end

