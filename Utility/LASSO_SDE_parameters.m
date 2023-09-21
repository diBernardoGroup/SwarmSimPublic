function [mu, theta, sigma] = LASSO_SDE_parameters(data, deltaT)
%LASSO_SDE_parameters Lasso estimation of the parameters of a SDE in the form
%         dX = theta * (mu - X) * dt + sigma * dW
%     where dW is gaussian white noise.
%     
%     The algorithm is described in "Calibrating the Ornstein-Uhlenbeck (Vasicek) model" by van den Berg (2011).

number_of_series = size(data,2);
mu=nan(1,number_of_series);
theta=nan(1,number_of_series);
sigma=nan(1,number_of_series);

for i=1:size(data,2)
    x = data(1:end-1,i);
    y = data(2:end,i);
    assert(length(x)==length(y))
    n=length(x);
    
    % leas sqaure regression
    X = [x];
    [B,FitInfo] = lasso(X,y,'CV',10 ,'Standardize', false);
    p=B(:,FitInfo.IndexMinMSE);
    a = p(1);
    b = FitInfo.Intercept(FitInfo.IndexMinMSE);
    residuals = y - (a*x + b);
    
    % compute estimated parameters
    mu(i) = b/(1-a);
    theta(i) = -1/deltaT * log(a);
    sigma(i) = std(residuals) * sqrt(-2*log(a) / (1-a^2) / deltaT);
end

end

