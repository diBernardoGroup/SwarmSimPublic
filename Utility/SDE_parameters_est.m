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
if isa(u,'cell')
    number_of_inputs = size(u{1},2);
else
    number_of_inputs = size(u,2);
end

mu=nan(number_of_series,1);
theta=nan(number_of_series,1);
sigma=nan(number_of_series,1);
alpha=nan(number_of_series,number_of_inputs);


% epsilon = 0.001;
figure
series_to_plot = round(linspace(1,number_of_series,min(number_of_series,5)));
ii=1;
for i=1:number_of_series
    nan_ids = isnan(x(:,i));
    d = x(~nan_ids,i);
    if isa(u,'cell')
        u_current = u{i}(~nan_ids,:);
        u_is_empty = all(u{i}==0,'all');
    else
        u_current = u(~nan_ids,:);
        u_is_empty = all(u==0,'all');
    end
    x_old = d(1:end-1);
    x_new = d(2:end);
    u_current = u_current(1:end-1,:);
    assert(length(x_old)==length(x_new))
    
    % estimate discrete time parameters
    if length(x_new) > 10 && (u_is_empty || all(var(u_current) > 0))
        % LASSO regression
        if strcmp(method,'LASSO')
            predictors = [x_old, u_current];
            [B,FitInfo] = lasso(predictors, x_new,'CV',10 ,'Standardize', false, 'Intercept',true);
            p=B(:,FitInfo.IndexMinMSE);
            a = p(1);
            b = p(2:end);
            c = FitInfo.Intercept(FitInfo.IndexMinMSE);
            
        % OLS regression
        elseif strcmp(method,'OLS')
            predictors = [x_old, u_current, ones(length(x_old),1)];
            p = predictors\x_new;
            a = p(1);
            b = p(2:end-1);
            c = p(end);
        
        % MATLAB DT grey box estimation
        elseif strcmp(method,'GreyBox')
            data = iddata(x_old, u_current, deltaT);
            sys = idnlgrey('discretizedSys',[1 number_of_inputs 1], {0, [0;0], 0}, x_old(1), deltaT);
            sys.Parameters(1).Name = 'a'; sys.Parameters(1).Minimum = 0; sys.Parameters(1).Maximum = 1;
            sys.Parameters(2).Name = 'b'; sys.Parameters(3).Name = 'c';
            %y = sim(sys, data);
            opt = nlgreyestOptions;
            %opt.Display = 'on';
            opt.EstimateCovariance = false;
            %opt.SearchOptions.FunctionTolerance = 10e-9;
            sys_id = nlgreyest(data, sys, opt);
            a = sys_id.Parameters(1).Value;
            b = sys_id.Parameters(2).Value;
            c = sys_id.Parameters(3).Value;
            
        % MATLAB direct CT grey box estimation
        elseif strcmp(method,'GreyBoxCT')
            data = iddata(x_old, u_current, deltaT);
            sys = idnlgrey('continouosSys',[1 number_of_inputs 1], {0, [0;0], 0}, x_old(1));
            sys.Parameters(1).Name = 'theta';   sys.Parameters(1).Minimum = 0; %sys.Parameters(1).Maximum = 1;
            sys.Parameters(2).Name = 'alpha'; 
            sys.Parameters(3).Name = 'mu';      %sys.Parameters(3).Minimum = 0;
            %y = sim(sys, data);
            opt = nlgreyestOptions;
            opt.EstimateCovariance = false;
            %opt.SearchOptions.FunctionTolerance = 10e-9;
            %opt.Display = 'on';
            sys_id = nlgreyest(data, sys, opt);
            t = sys_id.Parameters(1).Value;
            a = sys_id.Parameters(2).Value;
            m = sys_id.Parameters(3).Value;
            
            residuals = (x_new-x_old) - (t * (m - x_old) + u_current*a) * deltaT;
            s = std(residuals)/sqrt(deltaT);
            
            % plot
            if i >= series_to_plot(ii)
                subplot(1,length(series_to_plot),ii)
                hold on
                scatter(x_old, (x_new-x_old)/deltaT - u_current*a)
                plot(x_old, (t * (m-x_old)) )
                %axis([min(x,[],'all')*0.9, max(x,[],'all')*1.1, min(x,[],'all')*0.9, max(x,[],'all')*1.1])
                ii=ii+1;
            end
            
            theta(i) = t;
            alpha(i,:) = a;
            mu(i)    = m;
            sigma(i) = s;
        
        else
            error('Use a valid identification method!')
        end
        
        if exist('b','var')
            %compute continuous time parameters
            residuals = x_new - (a*x_old + u_current*b + c);
            mu(i) = real(c/(1-a));
            alpha(i,:) = real(inv(a-1)*log(a)/deltaT*b);
            if a ~= 0
                theta(i) = real(-1/deltaT * log(a));
                sigma(i) = real(std(residuals) * sqrt(-2*log(a) / (1-a^2) / deltaT));
            else
                theta(i) = 0;
                sigma(i) = 0;
            end
            
            % plot
            if i >= series_to_plot(ii)
                subplot(1,length(series_to_plot),ii)
                hold on
                scatter(x_old, x_new - u_current*b)
                plot(x_old, a*x_old + c)
                axis([min(x,[],'all')*0.9, max(x,[],'all')*1.1, min(x,[],'all')*0.9, max(x,[],'all')*1.1])
                ii=ii+1;
            end
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

