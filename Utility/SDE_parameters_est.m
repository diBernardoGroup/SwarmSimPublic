function [mu, theta, sigma, alpha] = SDE_parameters_est(x, u, deltaT, method, min_duration, no_mu)
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
    min_duration = 10 %[s]
    no_mu = false % fix mu to zero
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

% GreyBox estimation options
opt = nlgreyestOptions;
opt.EstimateCovariance = false;
%opt.Regularization.Lambda = 1;
%opt.SearchOptions.FunctionTolerance = 10e-9;
%opt.Display = 'on';

figure
series_to_plot = round(linspace(1,number_of_series,min(number_of_series,5)));
ii=1;
for i=1:number_of_series % for each time-series
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
    if length(x_new)*deltaT >= min_duration && (u_is_empty || all(var(u_current) > 0))
        % LASSO regression
        if strcmp(method,'LASSO')
            predictors = [x_old, u_current];
            [B,FitInfo] = lasso(predictors, x_new,'CV',10 ,'Standardize', true, 'Intercept',true);
            p=B(:,FitInfo.Index1SE);
            a = p(1);
            b = p(2:end);
            c = FitInfo.Intercept(FitInfo.Index1SE);
            %compute continuous time parameters
            residuals = x_new - (a*x_old + u_current*b + c);
            s=std(residuals);
            [theta(i), alpha(i,:), mu(i), sigma(i)] = dt2ct_parameters(a,b,c,s, deltaT);
            
            
        % OLS regression
        elseif strcmp(method,'OLS') || strcmp(method,'OLS+GB')
            if no_mu % fix mu to zero
                predictors = [x_old, u_current];
                p = predictors\x_new;
                a = p(1);
                b = p(2:end);
                c = 0;
            else
                predictors = [x_old, u_current, ones(length(x_old),1)];
                p = predictors\x_new;
                a = p(1);
                b = p(2:end-1);
                c = p(end);
            end
            %compute continuous time parameters
            residuals = x_new - (a*x_old + u_current*b + c);
            s=std(residuals);
            [theta(i), alpha(i,:), mu(i), sigma(i)] = dt2ct_parameters(a,b,c,s, deltaT);
            
            
        % MATLAB DT grey box estimation
        elseif strcmp(method,'GreyBoxDT')
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
            %compute continuous time parameters
            residuals = x_new - (a*x_old + u_current*b + c);
            s=std(residuals);
            [theta(i), alpha(i,:), mu(i), sigma(i)] = dt2ct_parameters(a,b,c,s, deltaT);
        end
        
        % MATLAB CT grey box estimation
        if strcmp(method,'GreyBoxCT') || strcmp(method,'OLS+GB')
            data = iddata(x_old, u_current, deltaT);
            sys = idnlgrey('continouosSys',[1 number_of_inputs 1], {1, [0;0], mean(x_old)}, x_old(1), 0);
            sys.Parameters(1).Name = 'theta';   sys.Parameters(1).Minimum = 0; %sys.Parameters(1).Maximum = 1;
            sys.Parameters(2).Name = 'alpha';
            sys.Parameters(3).Name = 'mu';      %sys.Parameters(3).Minimum = 0;
            if strcmp(method,'OLS+GB')
                sys.Parameters(1).Value = max(sys.Parameters(1).Minimum, theta(i));
                sys.Parameters(2).Value = alpha(i,:);
                %sys.Parameters(2).Value(1) = 0; sys.Parameters(2).Fixed(1) = true; % fix first component of alpha to zero
                sys.Parameters(3).Value = mu(i);
            end
            if no_mu % fix mu to zero
                sys.Parameters(3).Value = 0;
                sys.Parameters(3).Fixed = true; 
            end
            %y = sim(sys, data);
            opt.SearchOptions.FunctionTolerance = 10e-6;
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
        end
        
%         if exist('b','var')
%             % plot
%             if i >= series_to_plot(ii)
%                 subplot(1,length(series_to_plot),ii)
%                 hold on
%                 scatter(x_old, x_new - u_current*b)
%                 plot(x_old, a*x_old + c)
%                 axis([min(x,[],'all')*0.9, max(x,[],'all')*1.1, min(x,[],'all')*0.9, max(x,[],'all')*1.1])
%                 ii=ii+1;
%             end
%         end
    end
    
end

title(method)

end

