
theta=1;
mu=5;
sigma=0.5;

Tmax    = 100;
dT      = 0.01; % integration step
deltaT  = 1; % sampling time
Ntimes  = 10;

sde=hwv(theta,mu,sigma, 'StartState', 1); 
[paths,times]=simByEuler(sde, Tmax/dT, 'nTrials', Ntimes, 'DeltaTime', dT); 
%[paths,times]=simulate(sde, Tmax/dT, 'nTrials', Ntimes, 'DeltaTime', dT); 
data=squeeze(paths);

data=data(1:deltaT/dT:end,:);
times=times(1:deltaT/dT:end);

figure; plot(times, data)

[m_LASSO, t_LASSO, s_LASSO] = LASSO_SDE_parameters(data, deltaT);
[m_OLS, t_OLS, s_OLS]       = LS_SDE_parameters(data, deltaT);
[m_MLE, t_MLE, s_MLE]       = MLE_SDE_parameters(data, deltaT);

fprintf('Parameters: Tmax=%.1f\tdT=%.3f\tdeltaT=%.3f\tNtimes=%d\n', Tmax, dT, deltaT, Ntimes)
fprintf('Parameters: mean(data)=%.3f\tstd(data)=%.3f\n\n', mean(data(:)), std(data(:)))
fprintf('\ttheta\tmu\tsigma\n')
fprintf('TRUE\t%.2f\t%.2f\t%.2f\n',theta,mu,sigma)
fprintf('LASSO\t%.2f\t%.2f\t%.2f\n',mean(t_LASSO),  mean(m_LASSO),  mean(s_LASSO))
fprintf('OLS\t%.2f\t%.2f\t%.2f\n',  mean(t_OLS),    mean(m_OLS),    mean(s_OLS))
fprintf('MLE\t%.2f\t%.2f\t%.2f\n',  mean(t_MLE),    mean(m_MLE),    mean(s_MLE))

