function error = mape(F,A,option)
%MAPE Summary of this function goes here
%   Detailed explanation goes here

arguments
    F
    A
    option = 'default'
end

if strcmp(option,'default')
    error = mean(abs((A-F)./A));
elseif strcmp(option,'wMAPE')
    error = mean(abs(A-F))/mean(A);
elseif strcmp(option,'symmetric')
    error = 2*mean(abs(A-F)./(abs(A)+abs(F)));
end

end

