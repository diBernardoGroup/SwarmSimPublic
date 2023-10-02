function [mu, theta, sigma] = MLE_SDE_parameters(data, deltaT)
%MLE_SDE_PARAMETERS Maximum Likelihood Estimation of the parameters of a SDE in the form
%         dX = theta * (mu - X) * dt + sigma * dW
%     where dW is gaussian white noise.
%     
%     The algorithm is described in "Calibrating the Ornstein-Uhlenbeck (Vasicek) model" by van den Berg (2011).
    
number_of_series = size(data,2);
mu=nan(1,number_of_series);
theta=nan(1,number_of_series);
sigma=nan(1,number_of_series);

for i=1:number_of_series
    n=length(data)-1;
    x = data(1:end-1,i);
    y = data(2:end,i);
    assert(length(x)==length(y))
    
    % compute moments
    s_x = sum(x);
    s_y = sum(y);
    s_xx = sum(x.^2);
    s_yy = sum(y.^2);
    s_xy = sum(x .* y);
    
    % compute estimated parameters
    mu(i) = (s_y*s_xx - s_x*s_xy) / (n *(s_xx - s_xy) - (s_x^2 - s_x*s_y));
    theta(i) = -1/deltaT * log((s_xy - mu(i)*s_x - mu(i)*s_y + n*mu(i)^2)/(s_xx - 2*mu(i)*s_x + n*mu(i)^2));
    
    a = exp(-theta(i)*deltaT);
    sigma(i) = sqrt(2*theta(i)/(1-a^2) * 1/n * (s_yy - 2*a*s_xy + a^2*s_xx - 2*mu(i)*(1-a)*(s_y - a*s_x) + n*mu(i)^2*(1-a)^2));
end
end

