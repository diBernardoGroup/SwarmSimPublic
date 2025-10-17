function [mu, theta, sigma] = ARMAX_estimation(data, deltaT)

number_of_series = size(data,2);
mu=nan(1,number_of_series);
theta=nan(1,number_of_series);
sigma=nan(1,number_of_series);

for i=1:size(data,2)
    x = data(1:end-1,i);
    
    % leas sqaure regression
    IOdata = iddata(x,ones(size(x)),deltaT);
    sys = armax(IOdata, [1, 1, 1, 1],'Ts',deltaT);
    
    % compute estimated parameters
    mu(i) = b/(1-a);
    theta(i) = -1/deltaT * log(a);
    sigma(i) = std(residuals) * sqrt(-2*log(a) / (1-a^2) / deltaT);
end

end

