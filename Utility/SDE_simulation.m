clear


theta=1;
mu=5;
sigma=0.5;
alpha = 0.25;

Tmax    = 100;
dT      = 0.01; % integration step
deltaT  = 1; % sampling time
Ntimes  = 10;
times = [0:dT:Tmax];
%u = zeros(size(times));
%u = times>Tmax/3 & times<Tmax*2/3;
u = sin(times/10);

% sde=hwv(theta,mu,sigma, 'StartState', 1); 
% [paths,times]=simByEuler(sde, Tmax/dT, 'nTrials', Ntimes, 'DeltaTime', dT); 
% [paths,times]=simulate(sde, Tmax/dT, 'nTrials', Ntimes, 'DeltaTime', dT); 
% data=squeeze(paths);

X=nan(length(times),Ntimes);
X(1,:)=zeros(1,Ntimes);
for i=1:length(times)
    X(i+1,:)= X(i,:) + (theta*(mu-X(i,:)) + alpha*u(i))*dT + sigma* sqrt(dT) * randn(1,Ntimes);
end

X_data=X(1:deltaT/dT:end,:);
u_data=u(1:deltaT/dT:end)';
times=times(1:deltaT/dT:end);

figure; 
subplot(2,1,1)
plot(times, X_data)
subplot(2,1,2)
plot(times, u_data)

[m_LASSO, t_LASSO, s_LASSO, a_LASSO] = LASSO_SDE_parameters(X_data, u_data, deltaT);
[m_OLS, t_OLS, s_OLS]       = LS_SDE_parameters(X_data, deltaT);
[m_MLE, t_MLE, s_MLE]       = MLE_SDE_parameters(X_data, deltaT);

fprintf('Parameters: Tmax=%.1f\tdT=%.3f\tdeltaT=%.3f\tNtimes=%d\n', Tmax, dT, deltaT, Ntimes)
fprintf('Parameters: mean(data)=%.3f\tstd(data)=%.3f\n\n', mean(X_data(:)), std(X_data(:)))
fprintf('\ttheta\tmu\tsigma\talpha\n')
fprintf('TRUE\t%.2f\t%.2f\t%.2f\t%.2f\n',theta,mu,sigma,alpha)
fprintf('LASSO\t%.2f\t%.2f\t%.2f\t%.2f\n',mean(t_LASSO),  mean(m_LASSO),  mean(s_LASSO),  mean(a_LASSO))
fprintf('OLS\t%.2f\t%.2f\t%.2f\n',  mean(t_OLS),    mean(m_OLS),    mean(s_OLS))
fprintf('MLE\t%.2f\t%.2f\t%.2f\n',  mean(t_MLE),    mean(m_MLE),    mean(s_MLE))

