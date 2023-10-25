function [mu, theta, sigma] = LS_SDE_parameters(data, deltaT)
%LS_SDE_parameters Least Sqare estimation of the parameters of a SDE in the form
%         dX = theta * (mu - X) * dt + sigma * dW
%     where dW is gaussian white noise.
%     
%     The algorithm is described in "Calibrating the Ornstein-Uhlenbeck (Vasicek) model" by van den Berg (2011).

number_of_series = size(data,2);
mu=nan(1,number_of_series);
theta=nan(1,number_of_series);
sigma=nan(1,number_of_series);

for i=1:number_of_series
    nan_ids = isnan(data(:,i));
    d = data(~nan_ids,i);
    
    x = d(1:end-1);
    y = d(2:end);
    assert(length(x)==length(y))
    
    if length(x) > 2
        % leas sqaure regression
        p = polyfit(x,y,1);
        a = p(1);
        b = p(2);
        residuals = y - (a*x + b);

        % compute estimated parameters
        mu(i) = real(b/(1-a));
        theta(i) = real(-1/deltaT * log(a));
        sigma(i) = real(std(residuals) * sqrt(-2*log(a) / (1-a^2) / deltaT));
    end
end

end

