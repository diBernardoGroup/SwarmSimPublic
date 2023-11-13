function [mu, theta, sigma, alpha] = SDE_parameters_est(x, u, deltaT, method)
%SDE_parameters_est Regularized estimation of the parameters of a SDE in the form
%         dX = [ theta * (mu - X) + alpha * u] * dt + sigma * dW
%     where dW is gaussian white noise.
%     
%     The algorithm is described in "Calibrating the Ornstein-Uhlenbeck (Vasicek) model" by van den Berg (2011).

arguments
    x
    u
    deltaT = 1
    method = 'OLS'
end

number_of_series = size(x,2);
mu=nan(number_of_series,1);
theta=nan(number_of_series,1);
sigma=nan(number_of_series,1);
alpha=nan(number_of_series,size(u,2));

% epsilon = 0.001;
figure
series_to_plot = round(linspace(1,number_of_series,min(number_of_series,5)));
ii=1;
for i=1:number_of_series
    nan_ids = isnan(x(:,i));
    d = x(~nan_ids,i);
    u_current = u(~nan_ids,:);

    x_old = d(1:end-1);
    x_new = d(2:end);
    u_current = u_current(1:end-1,:);
    assert(length(x_old)==length(x_new))
    
    % LASSO regression
    % invert sampling equation
    if length(x_new) > 10
        if strcmp(method,'LASSO')
            predictors = [x_old, u_current];
            [B,FitInfo] = lasso(predictors, x_new,'CV',10 ,'Standardize', false, 'Intercept',true);
            p=B(:,FitInfo.IndexMinMSE);
            a = p(1);
            c = p(2:end);
            b = FitInfo.Intercept(FitInfo.IndexMinMSE);
        elseif strcmp(method,'OLS')
            predictors = [ones(length(x_old),1), x_old, u_current];
            p = predictors\x_new;
            a = p(2);
            b = p(1);
            c = p(3:end);
        end

        residuals = x_new - (b + a*x_old + u_current*c);
%         residuals_constant = x_new - mean(x_new);    
%         if norm(residuals_constant) < norm(residuals) + 0.01
%             residuals = residuals_constant;
%             b=mean(x_new);
%             a=0;
%         end
        
        if i >= series_to_plot(ii)
            subplot(1,length(series_to_plot),ii)
            hold on
            scatter(x_old, x_new - u_current*c)
            plot(x_old, b + a*x_old)
            axis([min(x,[],'all')*0.9, max(x,[],'all')*1.1, min(x,[],'all')*0.9, max(x,[],'all')*1.1])
            ii=ii+1;
        end
        
        %compute continuous time parameters
        mu(i) = real(b/(1-a));
        alpha(i,:) = real(inv(a-1)*log(a)/deltaT*c);
        if a ~= 0
            theta(i) = real(-1/deltaT * log(a));
            sigma(i) = real(std(residuals) * sqrt(-2*log(a) / (1-a^2) / deltaT));
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

title(method)

end

