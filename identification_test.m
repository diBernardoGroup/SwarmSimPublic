%% new identification test

theta=1;
mu=0;
sigma=0.5;

Tmax    = 100;
dT      = 0.01; % integration step
deltaT  = 0.5; % sampling time
Ntimes  = 1;

sde=hwv(theta,mu,sigma, 'StartState', 1); 
[paths,times]=simByEuler(sde, Tmax/dT, 'nTrials', Ntimes, 'DeltaTime', dT); 
%[paths,times]=simulate(sde, Tmax/dT, 'nTrials', Ntimes, 'DeltaTime', dT); 
data=squeeze(paths);

data=data(1:deltaT/dT:end,:);
times=times(1:deltaT/dT:end);

figure; plot(times, data)

IOdata = iddata(data,zeros(length(times),1),deltaT)

% define structure of the model
order=1;
A=ones(order);
B=ones(order,1);
C=zeros(1,order); C(1)=1;
identified_sys = idss(A,B,C,0,ones(order,1),zeros(order,1),deltaT)
identified_sys.Structure.C.Free=false;

% identify parameters
identified_sys = ssest(IOdata, identified_sys)

sys_c = d2c(identified_sys)


