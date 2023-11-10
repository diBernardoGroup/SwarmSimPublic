function [mu, theta, sigma] = LS_SDE_parameters(x, deltaT)
%LS_SDE_parameters Least Sqare estimation of the parameters of a SDE in the form
%         dX = theta * (mu - X) * dt + sigma * dW
%     where dW is gaussian white noise.
%     
%     The algorithm is described in "Calibrating the Ornstein-Uhlenbeck (Vasicek) model" by T. van den Berg (2011).

number_of_series = size(x,2);
mu=nan(1,number_of_series);
theta=nan(1,number_of_series);
sigma=nan(1,number_of_series);

for i=1:number_of_series
    nan_ids = isnan(x(:,i));
    d = x(~nan_ids,i);
    
    x_old = d(1:end-1);
    x_new = d(2:end);
    assert(length(x_old)==length(x_new))
    
    if length(x_new) > 2
        % leas sqaure regression
        p = polyfit(x_old,x_new,1);
        a = p(1);
        b = p(2);
        residuals = x_new - (a*x_old + b);

        % compute estimated parameters
        mu(i) = real(b/(1-a));
        theta(i) = real(-1/deltaT * log(a));
        sigma(i) = real(std(residuals) * sqrt(-2*log(a) / (1-a^2) / deltaT));
    end
end

end

