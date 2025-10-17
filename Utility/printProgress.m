function printProgress(current,total)
%PRINTPROGERESS Summary of this function goes here
%   Detailed explanation goes here
    fprintf(repmat('\b', 1, 7*(current>1)));
    fprintf('%6.1f%%',current/total*100)
    fprintf(repmat('\n', 1, (current==total)));
end

