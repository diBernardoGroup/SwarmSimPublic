function [mu, theta, sigma, alpha] = SDE_parameters_est(x, u, deltaT, method, min_duration, no_mu, model, limits, make_plots)
%SDE_parameters_est Regularized estimation of the parameters of a SDE in the form
%         dX = [ theta * (mu - X) + alpha * u] * dt + sigma * dW
%     where dW is gaussian white noise.
%
%     The algorithm is described in "Calibrating the Ornstein-Uhlenbeck (Vasicek) model" by van den Berg (2011).

arguments
    x
    u
    deltaT = 1
    method (1,:) char {mustBeMember(method,{'LASSO','OLS','OLS+GB','GreyBoxDT','GBDT','GreyBoxCT','GBCT','GB','GA'})} ='OLS'
    min_duration = 10                   %[s]
    no_mu = false                       % fix mu to zero
    model = ''
    limits= [[0;inf], inf(2,3).*[-1;1]]
    make_plots = false
end

number_of_series = size(x,2);
if isempty(limits)
    limits = [[0;inf], inf(2,3).*[-1;1]];
end
if no_mu
    limits(:,4) = [0;0]; 
end
if isa(u,'cell')
    number_of_inputs = size(u{1},2);
    u_is_empty = all([u{:}]==0 | isnan([u{:}]),'all');
else
    number_of_inputs = size(u,2);
    u_is_empty = all(u==0,'all');           % no inputs
end

mu=nan(number_of_series,1);
theta=nan(number_of_series,1);
sigma=nan(number_of_series,1);
alpha=nan(number_of_series,number_of_inputs);

% GreyBox estimation options
optGB = nlgreyestOptions;
optGB.EstimateCovariance = false;
optGB.SearchOptions.FunctionTolerance = 10e-6;

% Genetic Algorithm options
optGA = optimoptions('ga','PopulationSize',200,'MaxStallGenerations',100,'MaxGenerations',500,'display','off');

if make_plots; figure; end
series_to_plot = round(linspace(1,number_of_series,min(number_of_series,5)));
ii=1;
fprintf('> %s identification',method)

for i=1:number_of_series % for each time-series
    printProgress(i,number_of_series)
    
    nan_ids = isnan(x(:,i));
    d = x(~nan_ids,i);                      % d is the trajectory of x (when it exist)
    
    if isa(u,'cell')
        u_current = u{i}(~nan_ids,:);       % u_current is the control input
    else
        u_current = u(~nan_ids,:);          % u_current is the control input
    end
    
    u_current = u_current(1:end-1,:);       % Current input
    x_old = d(1:end-1);                     % Current state
    x_new = d(2:end);                       % Next state
    assert(length(x_old)==length(x_new) && length(x_new)==size(u_current,1))
    
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
        elseif strcmp(method,'GreyBoxDT') || strcmp(method,'GBDT')
            data = iddata(x_old, u_current, deltaT);
            sys = idnlgrey('discretizedSys',[1 number_of_inputs 1], {0, [0;0], 0}, x_old(1), deltaT);
            sys.Parameters(1).Name = 'a'; sys.Parameters(1).Minimum = 0; sys.Parameters(1).Maximum = 1;
            sys.Parameters(2).Name = 'b'; 
            sys.Parameters(3).Name = 'c';
            sys_id = nlgreyest(data, sys, optGB);
            a = sys_id.Parameters(1).Value;
            b = sys_id.Parameters(2).Value;
            c = sys_id.Parameters(3).Value;
            %compute continuous time parameters
            residuals = x_new - (a*x_old + u_current*b + c);
            s=std(residuals);
            [theta(i), alpha(i,:), mu(i), sigma(i)] = dt2ct_parameters(a,b,c,s, deltaT);
           
        % Genetic Algorithm
        elseif strcmp(method,'GA')
            pars = ga(@(a) model(a, d, u_current, deltaT),4,[],[],[],[],limits(1,:),limits(2,:),[],optGA);
            %pars = ga(@(a) id_cost(a, d, u_current, deltaT, model),4,[],[],[],[],limits(1,:),limits(2,:),[],optGA);

            theta(i) = pars(1);
            alpha(i,1) = pars(2);
            alpha(i,2) = pars(3);
            mu(i) = pars(4);
            
            % comute residuals to estimate sigma (noise amplitude)
            residuals = (x_new-x_old) - (pars(1) * (pars(4) - x_old) + sign(x_old).*u_current*[pars(2);pars(3)]) * deltaT;
            s = std(residuals)/sqrt(deltaT);
            sigma(i) = s ;        
        end
        
        % MATLAB CT grey box estimation
        if strcmp(method,'GreyBoxCT') || strcmp(method,'GBCT') || strcmp(method,'GB') || strcmp(method,'OLS+GB')
            data = iddata(x_old, u_current, deltaT);
            sys = idnlgrey(model,[1 number_of_inputs 1], {1, [0;0], mean(x_old)}, x_old(1), 0);
            sys.Parameters(1).Name = 'theta';   sys.Parameters(1).Minimum = 0;
            sys.Parameters(2).Name = 'alpha';
            sys.Parameters(3).Name = 'mu';
            if strcmp(method,'OLS+GB')
                sys.Parameters(1).Value = max(sys.Parameters(1).Minimum, theta(i));
                sys.Parameters(2).Value = alpha(i,:);
                sys.Parameters(3).Value = mu(i);
            end
            l=1;
            for p=1:min(size(limits,2), length(sys.Parameters))
                lim_l = limits(1,l:l+length(sys.Parameters(p).Minimum)-1)';
                lim_u = limits(2,l:l+length(sys.Parameters(p).Minimum)-1)';
                sys.Parameters(p).Value = max(sys.Parameters(p).Value, lim_l);
                sys.Parameters(p).Value = min(sys.Parameters(p).Value, lim_u);
                
                sys.Parameters(p).Minimum = lim_l;
                sys.Parameters(p).Maximum = lim_u;
                
                sys.Parameters(p).Fixed = sys.Parameters(p).Maximum == sys.Parameters(p).Minimum;
                
                l=l+length(sys.Parameters(p).Maximum);
            end
            sys_id = nlgreyest(data, sys, optGB);
            t = sys_id.Parameters(1).Value;
            a = sys_id.Parameters(2).Value;
            m = sys_id.Parameters(3).Value;
            
            residuals = (x_new-x_old) - (t * (m - x_old) + u_current*a) * deltaT;
            s = std(residuals)/sqrt(deltaT);
            
            theta(i) = t;
            alpha(i,:) = a;
            mu(i)    = m;
            sigma(i) = s;
        end
        
        % plot
        if make_plots && i >= series_to_plot(ii)
            subplot(1,length(series_to_plot),ii)
            hold on
            scatter(x_old, (x_new-x_old)/deltaT - u_current*alpha(i,:)')
            plot(x_old, (theta(i) * (mu(i)-x_old)) )
            ii=ii+1;
        end
        
        
    end
    
    if make_plots; title(method); end
    
end

