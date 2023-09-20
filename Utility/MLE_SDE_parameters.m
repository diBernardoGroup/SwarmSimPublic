function [mu, theta, sigma] = MLE_SDE_parameters(data, deltaT)
%MLE_SDE_PARAMETERS Estimate parameters of a SDE in the form
%         dX = theta * (mu - X) * dt + sigma * dW
%     where dW is gaussian white noise.
%     
%     The algorithm is described in "Calibrating the Ornstein-Uhlenbeck (Vasicek) model" by van den Berg (2011).
    
    n=length(data);
    x = data(1:end-1);
    y = data(2:end);
    assert(length(x)==length(y))
    
    % compute moments
    s_x = sum(x);
    s_y = sum(y);
    s_xx = sum(x.^2);
    s_yy = sum(y.^2);
    s_xy = sum(x .* y);
    
    % compute estimated parameters
    mu = (s_y*s_xx - s_x*s_xy) / (n *(s_xx - s_xy) - (s_x^2 - s_x*s_y));
    theta = -1/deltaT * log((s_xy - mu*s_x - mu*s_y + n*mu^2)/(s_xx - 2*mu*s_x + n*mu^2));
    
    a = exp(-theta*deltaT);
    sigma = sqrt(2*theta/(1-a^2) * 1/n * (s_yy - 2*a*s_xy + a^2*s_xx - 2*mu*(1-a)*(s_y - a*s_x) + n*mu.^2*(1-a)^2));

end

