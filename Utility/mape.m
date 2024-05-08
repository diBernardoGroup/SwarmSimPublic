function error = mape(x,x_ref,option)
%MAPE Summary of this function goes here
%   Detailed explanation goes here

arguments
    x
    x_ref
    option = 'default'
end

if strcmp(option,'default')
    error = mean(abs((x_ref-x)./x_ref));
elseif strcmp(option,'wMAPE')
    error = mean(abs(x_ref-x))/mean(x_ref);
elseif strcmp(option,'symmetric')
    error = 2*mean(abs(x_ref-x)./(abs(x_ref)+abs(x)));
end

end

