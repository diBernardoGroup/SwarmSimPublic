function [mu, theta, sigma, alpha] = LASSO_SDE_parameters(x, u, deltaT)
%LASSO_SDE_parameters Lasso estimation of the parameters of a SDE in the form
%         dX = [ theta * (mu - X) + alpha * u] * dt + sigma * dW
%     where dW is gaussian white noise.
%     
%     The algorithm is described in "Calibrating the Ornstein-Uhlenbeck (Vasicek) model" by van den Berg (2011).

number_of_series = size(x,2);
mu=nan(1,number_of_series);
theta=nan(1,number_of_series);
sigma=nan(1,number_of_series);
alpha=nan(1,number_of_series);

if length(u)==size(x,1)
   u=u(1:end-1); 
end

epsilon = 0.001;

for i=1:size(x,2)
    x_old = x(1:end-1,i);
    x_new = x(2:end,i);
    assert(length(x_old)==length(x_new))
    
    % LASSO regression
    % invert sampling equation
    X = [x_old, u];
    [B,FitInfo] = lasso(X, x_new,'CV',10 ,'Standardize', false, 'Intercept',true);
    p=B(:,FitInfo.IndexMinMSE);
    a = p(1);
    c = p(2);
    b = FitInfo.Intercept(FitInfo.IndexMinMSE);
    
    residuals = x_new - (b + a*x_old + c*u);
    residuals_constant = x_new - mean(x_new);    
    if norm(residuals_constant) < norm(residuals) + 0.01
        residuals = residuals_constant;
        b=mean(x_new);
        a=0;
    end
    
    %compute estimated parameters
    mu(i) = b/(1-a);
    alpha(i) = inv(a-1)*log(a)/deltaT*c;
    if a ~= 0
        theta(i) = -1/deltaT * log(a);
        sigma(i) = std(residuals) * sqrt(-2*log(a) / (1-a^2) / deltaT);
    else
        theta(i) = 0;
        sigma(i) = 0;
    end
    
%     x_old = x(1:end-1,i);
%     x_new = x(2:end,i) - x;
%     assert(length(x_old)==length(x_new))
%     
%     % invert forward Euler integration equation
%     X = [x];
%     [B,FitInfo] = lasso(X,x_new,'CV',10 ,'Standardize', false, 'Intercept',true);
%     p=B(:,FitInfo.IndexMinMSE);
%     a = p(1);
%     b = FitInfo.Intercept(FitInfo.IndexMinMSE);
%     
%     residuals = x_new - (a*x + b);
%     residuals_constant = x_new - mean(x_new);    
% %     if norm(residuals_constant) < norm(residuals) + 0.01
% %         residuals = residuals_constant;
% %         b=mean(x_new);
% %         a=0;
% %     end
%     
%     % compute estimated parameters
%     theta(i) = -a/deltaT;
%     if abs(a) < epsilon
%         mu(i) = mean(x);
%     else
%         mu(i) = -b/a;
%     end
%     sigma(i) = std(residuals) / sqrt(deltaT);
    
end

end

