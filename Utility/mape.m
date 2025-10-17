function error = mape(x,x_ref,option)
%MAPE Summary of this function goes here
%   Detailed explanation goes here

arguments
    x
    x_ref
    option = 'default'
end

if strcmp(option,'default')     % MAPE
    error = mean(abs((x_ref-x)./x_ref),'omitnan');
elseif strcmp(option,'wMAPE')   % Weighted MAPE (on refernece signal average)
    error = mean(abs(x_ref-x),'omitnan')/mean(x_ref,'omitnan');
elseif strcmp(option,'symmetric') % Symmetric MAPE (on both signals averages)
    error = 2*mean(abs(x_ref-x)./(abs(x_ref)+abs(x)),'omitnan');
end

end

