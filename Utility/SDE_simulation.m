
theta=1;
mu=5;
sigma=0.5;

deltaT = 0.25;
Ntimes =10;

sde=hwv(theta,mu,sigma); 
[paths,times]=simulate(sde, 100/deltaT, 'nTrials', Ntimes, 'DeltaTime', deltaT); 
data=squeeze(paths);

plot(times, data)

[m_LASSO, t_LASSO, s_LASSO] = LASSO_SDE_parameters(data, deltaT);
[m_OLS, t_OLS, s_OLS]       = LS_SDE_parameters(data, deltaT);
[m_MLE, t_MLE, s_MLE]       = MLE_SDE_parameters(data, deltaT);

fprintf('\ttheta\tmu\tsigma\n')
fprintf('TRUE\t%.3f\t%.3f\t%.3f\n',theta,mu,sigma)
fprintf('LASSO\t%.3f\t%.3f\t%.3f\n',mean(t_LASSO),  mean(m_LASSO),  mean(s_LASSO))
fprintf('OLS\t%.3f\t%.3f\t%.3f\n',  mean(t_OLS),    mean(m_OLS),    mean(s_OLS))
fprintf('MLE\t%.3f\t%.3f\t%.3f\n',  mean(t_MLE),    mean(m_MLE),    mean(s_MLE))

