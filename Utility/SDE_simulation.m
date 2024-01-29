clear
close all


theta   = 0.5;
mu      = 5;
sigma   = 0.1;
alpha   = -0.3;
beta    = -3;

Tmax    = 100;
dT      = 0.01; % integration step
deltaT  = 1; % sampling time
Ntimes  = 10;
times = [0:dT:Tmax];
% u = zeros(size(times));            % no input
u = times>Tmax/3 & times<Tmax*2/3;  % step
% u = times*Tmax;                     % ramp
% u = sin(times/10);                  % sine wave
% u = mod(times,10)/10;               % sawtooth wave

u_dot = [diff(u),0]/dT;
X=nan(length(times),Ntimes);
X(1,:)=rand(1,Ntimes)*2*mu;
for i=1:length(times)-1
    X(i+1,:)= X(i,:) + (theta*(mu-X(i,:)) + alpha*u(i) + beta*u_dot(i) )*dT + sigma* sqrt(dT) * randn(1,Ntimes);
end

X_data=X(1:deltaT/dT:end,:);
u_data=u(1:deltaT/dT:end)';
u_dot_data=[diff(u_data);0]/deltaT;
t_data=times(1:deltaT/dT:end);

[m_LASSO, t_LASSO, s_LASSO, a_LASSO]= SDE_parameters_est(X_data, [u_data, u_dot_data], deltaT, 'LASSO');
[m_OLS, t_OLS, s_OLS, a_OLS]        = SDE_parameters_est(X_data, [u_data, u_dot_data], deltaT, 'OLS');
[m_Grey, t_Grey, s_Grey, a_Grey]    = SDE_parameters_est(X_data, [u_data, u_dot_data], deltaT, 'GreyBox');
[m_GreyCT, t_GreyCT, s_GreyCT, a_GreyCT] = SDE_parameters_est(X_data, [u_data, u_dot_data], deltaT, 'GreyBoxCT');
[m_MLE, t_MLE, s_MLE]               = MLE_SDE_parameters(X_data, deltaT);

fprintf('Parameters: Tmax=%.1f\tdT=%.3f\tdeltaT=%.3f\tNtimes=%d\n', Tmax, dT, deltaT, Ntimes)
fprintf('Parameters: mean(data)=%.3f\tstd(data)=%.3f\n\n', mean(X_data(:)), std(X_data(:)))
fprintf('\ttheta\tmu\tsigma\talpha\tbeta\n')
fprintf('TRUE\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',theta,mu,sigma,alpha,beta)
fprintf(join(['LASSO\t%.2f\t%.2f\t%.2f',repmat("%.2f",1,length(mean(a_LASSO))),'\n'],'\t'),   mean(t_LASSO),   mean(m_LASSO),   mean(s_LASSO), mean(a_LASSO))
fprintf(join(['OLS\t%.2f\t%.2f\t%.2f',repmat("%.2f",1,length(mean(a_OLS))),'\n'],'\t'),       mean(t_OLS),     mean(m_OLS),     mean(s_OLS),   mean(a_OLS))
fprintf(join(['GreyBox\t%.2f\t%.2f\t%.2f',repmat("%.2f",1,length(mean(a_Grey))),'\n'],'\t'),  mean(t_Grey),    mean(m_Grey),    mean(s_Grey),  mean(a_Grey))
fprintf(join(['GB CT\t%.2f\t%.2f\t%.2f',repmat("%.2f",1,length(mean(a_GreyCT))),'\n'],'\t'),  mean(t_GreyCT),  mean(m_GreyCT),  mean(s_GreyCT),mean(a_GreyCT))
%fprintf('MLE\t%.2f\t%.2f\t%.2f\n',  mean(t_MLE),    mean(m_MLE),    mean(s_MLE))

%% simulate average identified system
x_sim=nan(length(times),1);
x_sim_OLS(1)=mean(m_OLS);
x_sim_Grey(1)=mean(m_Grey);
x_sim_GB_CT(1)=mean(m_GreyCT);
for i=1:length(times)-1
    x_sim_OLS(i+1) = x_sim_OLS(i) + (mean(t_OLS)*(mean(m_OLS)-x_sim_OLS(i)) + mean(a_OLS(:,1))*u(i) + mean(a_OLS(:,2))*u_dot(i) )*dT;
    x_sim_Grey(i+1)= x_sim_Grey(i) + (mean(t_Grey)*(mean(m_Grey)-x_sim_Grey(i)) + mean(a_Grey(:,1))*u(i) + mean(a_Grey(:,2))*u_dot(i) )*dT;
    x_sim_GB_CT(i+1)= x_sim_GB_CT(i) + (mean(t_GreyCT)*(mean(m_GreyCT)-x_sim_GB_CT(i)) + mean(a_GreyCT(:,1))*u(i) + mean(a_GreyCT(:,2))*u_dot(i) )*dT;
end

%% PLOTS
figure; 
subplot(3,1,1)
hold on
plot(times, X, color=[0.5,0.5,0.5])
plot(times, x_sim_OLS, 'r', LineWidth=1)
plot(times, x_sim_Grey, 'g', LineWidth=1)
plot(times, x_sim_GB_CT, LineWidth=1)
subplot(3,1,2)
hold on
plot(times, u)
plot(t_data, u_data,'--')
subplot(3,1,3)
hold on
plot(times, u_dot)
plot(t_data, u_dot_data,'--')



