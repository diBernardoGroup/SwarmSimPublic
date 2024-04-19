function tracking = getLastTracking(experimentDir)
%GETLASTTRACKING Summary of this function goes here
%   Detailed explanation goes here
all_tracking=dir(fullfile(experimentDir,'tracking_20*'));

assert(~isempty(all_tracking), char(['No tracking found in ',char(experimentDir)]))

tracking=all_tracking(end).name;
end

