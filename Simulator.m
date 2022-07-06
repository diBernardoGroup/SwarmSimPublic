function [T_r, success, final_e_theta, final_e_L, finalGRadial, finalGNormal, stopTime] = Simulator(x0, LinkNumber, GainRadialDefault, GainNormalDefault, regularity_thresh, compactness_thresh, Tmax, sigma, drawON, getMetrics, IntFunctionStruct, MaxSensingRadius, alpha, beta, dynamicLattice, AgentsRemoval)
%
%Simulator Executes a complete simulation of the swarm.
%   This function is called by a launcher script (Launcher, BruteForceTuning, ...).
%
%   [T_r, success, final_e_theta, final_e_L, finalGRadial, finalGNormal, stopTime] =...
%   ... Simulator(x0, LinkNumber, GainRadialDefault, GainNormalDefault, regularity_thresh, compactness_thresh,...
%   ... Tmax, sigma, drawON, getMetrics, IntFunctionStruct, MaxSensingRadius, alpha, beta, dynamicLattice, AgentsRemoval)
%
%   Inputs:
%       x0 are the initial positions of the agents (Nx2 matrix)
%       LinkNumber is the desired number of links per agent (6=triangular
%           lattice, 4=square lattice, 3=hexagonal lattice) (scalar)
%       GainRadialDefault is the default value of G_radial (scalar)
%       GainNormalDefault is the default value of G_normal (scalar)
%       regularity_thresh is the threshold value for regularity metrics (e^*_theta) (scalar)
%       compactness_thresh is the threshold value for compactness metrics (e^*_L) (scalar)
%       Tmax is the maximum simulation time (scalar)
%       sigma is the standard deviation of noise (scalar)
%       drawON draw swarm during simulation (bool)
%       getMetrics acquire metrics during the simulation (getMetrics=false
%           discard settling times and stop times) (bool)
%       IntFunctionStruct description and parameters of the radial
%           interaction function (struct)
%       MaxSensingRadius is the sensing radius of the agents (R_s)
%       alpha and beta are the adaptation gains for the adaptive control (scalar)
%       dynamicLattice change lattice during the simulation (bool)
%       AgentsRemoval randomly remove agents during the simulation (bool)
%
%   Outputs:
%       T_r is the settling time (scalar)
%       success is true if the metrics are below the respective thresholda (bool)
%       final_e_theta final value of the regularity metrics (scalar)
%       final_e_L final value of the compactness metrics (scalar)
%       finalGRadial final value of the radial gain (scalar)
%       finalGNormal final value of the normal gain (scalar)
%       stopTime final simulation instant (scalar)
%
%   See also: Launcher, BruteForceTuning, SequentialLauncher
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%

%% Check input parameters
assert(LinkNumber==6 | LinkNumber==4, "LinkNumber must be equal to 4 (square lattice) or 6 (triangular lattice)")

%% Instantiate Simulation Window
Max = 10;   % amplitude of the simulation plane
Min = -Max;

if drawON
    figure;
    axis('equal',[Min Max Min Max])
    yticks([-10 -5 0 5 10])
    xticks([-10 -5 0 5 10])
    set(gca,'FontSize',14)
    hold on
end

vMax = 5; % maximum speed of the agents

RMax = 1.1;             % maximum distance of adjacent agents
RMin = 0.6;             % minimum distance of adjacent agents
SensingNumber = inf;    % max number of neighbours to interact with

DeadZoneThresh=regularity_thresh; % amplitude of the dead zone in gain adaptation law

%% simulation parameters
InteractionFactor = 1; % fraction of agents to interact with [0,1] 
%(if InteractionFactor<1 a subset of agents is randomly selected by each agent at each update step to interact with)
deltaT = 0.01;      % forward Euler integration step
deltaSample = 0.25; % time step for metrics acquisition
screenTimes=[0 Tmax/2-deltaT Tmax/2 Tmax]; % specify time instants to get simulation frames

% steady state detection
stopSamples=ceil(10/deltaSample);               % number of consecutive samples to detect steady state 
regularity_ss_thresh=regularity_thresh/10;      % steady state detection threshold for the regularity metrics
compactness_ss_thresh=compactness_thresh/10;    % steady state detection threshold for the compactness metrics

%% Inizialization
x=x0;
N=size(x,1);

%Set spin for simulations with Spears' algorithm
spin=ones(N, 1);  
if strcmp(IntFunctionStruct.function,'Spears') && (LinkNumber==4 || LinkNumber == 3)
    spin(1:2:N, 1)=0;
end

if AgentsRemoval
    savedIndex = 1;
    xSaved = cell(length(screenTimes),1);
end


%% Preallocate variables
count=0;                                    % sampling iteration
TSample = 0:deltaSample:Tmax;               % sampling time instants
e_L=nan(size(TSample,2),3);                 % compactness metrics of the swarm
e_theta=nan(size(TSample,2),3);             % regularity metrics of the swarm
G_normal_vec=nan(size(TSample,2),3);        % average, min, and max values of the normal gain
G_radial_vec=nan(size(TSample));            % average radial gain
vMean=nan(size(TSample,2),3);               % average speed
G_radial = GainRadialDefault* ones(N,1);    % radial gains (G_r,i)
G_normal = GainNormalDefault* ones(N,1);    % normal gains (G_n,i)


disp(['Simulating ',IntFunctionStruct.function,' N=',num2str(N),' LinkNumber=',num2str(LinkNumber)])

%% Run Simulation
stopCondition=false;
t=0;
stopTime=nan;

while t<=Tmax && ~stopCondition
    
    if dynamicLattice
        if t>Tmax/3 && t<2*Tmax/3
            LinkNumber = 6;
        else
            LinkNumber = 4;
        end
    end
    
    % First Order Dynamics
    [v, links, ~, G_radial, G_normal, ~]= VFcontroller(x, [], 0, G_radial, G_normal, zeros(N,1), 0, 0, min(RMax,MaxSensingRadius), RMin, SensingNumber, InteractionFactor, LinkNumber, deltaT, DeadZoneThresh, IntFunctionStruct, spin, MaxSensingRadius, alpha, beta);
    x = SingleIntegrator(x, v, deltaT, vMax, sigma);
    
    if t>=TSample(count+1)
        count= count+1;
             
        if getMetrics || t>Tmax-2
            linkError=abs((LinkNumber-links))/LinkNumber;
            e_L(count,:)=[mean(linkError), quantile(linkError, 0.1), quantile(linkError, 0.9)];
            e_theta(count,1)=getAvgAngularErr_Spears(x, LinkNumber, RMin, RMax)*LinkNumber/pi;
            G_normal_vec(count,:)=[mean(G_normal), quantile(G_normal, 0.1), quantile(G_normal, 0.9)];
            G_radial_vec(count)=mean(G_radial);
            vMean(count,:)=[mean(vecnorm(v,2,2)), min(vecnorm(v,2,2)), max(vecnorm(v,2,2))];
           
            % steady-state detection (with vibration exclusion)
            if(count>stopSamples)
                if (abs(e_L(count-stopSamples:count,1)-e_L(count,1))< compactness_ss_thresh...
                        | abs(e_L(count-stopSamples:count,1)-e_L(count-1,1))< compactness_ss_thresh)...
                        & (abs(e_theta(count-stopSamples:count,1)-e_theta(count,1))< regularity_ss_thresh...
                        | abs(e_theta(count-stopSamples:count,1)-e_theta(count-1,1))< regularity_ss_thresh)...
                        & abs(G_normal_vec(count-stopSamples:count,1)-G_normal_vec(count,1))< G_normal_vec(count,1)*0.03
                    stopCondition=true;
                    stopTime=t;
                end
            end
            
        end
        
        % plot swarm
        if drawON
            plotSwarm(x, [], t, RMin,RMax, true, spin);
        end
        
        %random kill (remove a percentage of randomly selected agents)
        if t>Tmax/2-deltaT && size(x,1)==N && AgentsRemoval
            [x, spin]=randomKill(x, spin, 0.3);
        end
    end
    
    if ismember(t,screenTimes) && AgentsRemoval
        xSaved(savedIndex) = {x};
        savedIndex = savedIndex + 1;
        %plotSwarm(x, [], t, MinSensingRadius,MaxSensingRadius, true, spin);
    end
    
    t=t+deltaT;
end

% compute settling times
[T_r, T_theta, T_L, success, ~, ~] = settlingTimes(e_theta(:,1), e_L(:,1), regularity_thresh, compactness_thresh, TSample);
final_e_theta = mean(e_theta(count-ceil(1/deltaSample):count,1));
final_e_L = mean(e_L(count-ceil(1/deltaSample):count,1));
finalGRadial = mean(G_radial_vec(count-ceil(1/deltaSample):count),1);
finalGNormal = mean(G_normal_vec(count-ceil(1/deltaSample):count,1));


%% PLOTS

if drawON
    plotSwarm(x,[],t-deltaT,RMin,RMax,false, spin);
     
    figure %METRICS
    subplot(2,1,1)
    title('Metrics','FontSize',16)
    hold on
    e_theta_line=plot(TSample, e_theta(:,1),'b','LineWidth',3);
    e_theta_star_line=yline(regularity_thresh, 'b--','LineWidth',2);
    if ~isnan(T_theta); T_theta_line=xline(T_theta,'b--', 'Interpreter','latex','FontSize',16); end
    if AgentsRemoval; xline(Tmax/2,'b--', 'Interpreter','latex','FontSize',16); end
    axis([-inf inf 0 0.6])
    yticks([0 0.1 regularity_thresh 0.3 0.4 0.5 0.6])
    legend([e_theta_line, e_theta_star_line],{'$\bar{e}_{\theta}$','$e_{\theta}^*$'}, 'Interpreter','latex','FontSize',22)
    legend boxoff
    set(gca,'FontSize',14)
    
    subplot(2,1,2)
    e_L_line=plot(TSample,e_L(:,1), 'r','LineWidth',3);
    e_L_star_line=yline(compactness_thresh,'r--','LineWidth',2);
    if ~isnan(T_L); T_L_line=xline(T_L,'r--','Interpreter','latex', 'FontSize',16); end
    if AgentsRemoval; xline(Tmax/2,'r--', 'Interpreter','latex','FontSize',16); end
    yticks([0 0.15 compactness_thresh 0.45 0.6 0.75 0.9])
    axis([-inf inf 0 0.9])
    legend([e_L_line, e_L_star_line],{'$\bar{e}_{L}$','$e_{L}^*$'},'Interpreter','latex','FontSize',22)
    legend boxoff
    xlabel('t','FontSize',16)
    set(gca,'FontSize',14)
    
    if alpha~=0 || beta~=0 %GAINS
    figure 
    GTline=plotWithShade(TSample,G_normal_vec(:,1),G_normal_vec(:,2),G_normal_vec(:,3), 'b', 0.1);
    
    hold on
    GNline=plot(TSample, G_radial_vec,'LineWidth',3);
    if AgentsRemoval; xline(Tmax/2,'k--', 'Interpreter','latex','FontSize',16); end
    title('Gains','FontSize',10)
    legend([GTline GNline], {'G_n','G_r'},'FontSize',16)
    xlabel('t','FontSize',16)
    axis([-inf inf 0 max(G_radial_vec)+1])
    set(gca,'FontSize',14)
    end
    
end

if AgentsRemoval
    figure
    for i=1:length(xSaved)
        subplot(1,length(xSaved),i)
        axis('equal',[Min Max Min Max])
        yticks([-10 -5 0 5 10])
        xticks([-10 -5 0 5 10])
        set(gca,'FontSize',14)
        hold on
        plotSwarm(cell2mat(xSaved(i)),[],screenTimes(i),RMin,RMax,false, ones(size(cell2mat(xSaved(i)),1),1));
    end
    figure
    for i=1:length(xSaved)
        subplot(2,2,i)
        axis('equal',[Min Max Min Max])
        yticks([-10 -5 0 5 10])
        xticks([-10 -5 0 5 10])
        set(gca,'FontSize',14)
        hold on
        plotSwarm(cell2mat(xSaved(i)),[],screenTimes(i),RMin,RMax,false, ones(size(cell2mat(xSaved(i)),1),1));
    end
    for i=1:length(xSaved)
        figure
        axis('equal',[Min Max Min Max])
        yticks([-10 -5 0 5 10])
        xticks([-10 -5 0 5 10])
        set(gca,'FontSize',14)
        hold on
        plotSwarm(cell2mat(xSaved(i)),[],screenTimes(i),RMin,RMax,false, ones(size(cell2mat(xSaved(i)),1),1));
        title('')
    end
end


if dynamicLattice
    figure
    subplot(3,1,1)
    hold on
    e_theta_line=plot(TSample, e_theta(:,1), 'b', 'LineWidth',2);
    e_theta_star_line=yline(regularity_thresh, 'b--','LineWidth',2);
    legend('$\bar{e}_{\theta}$', '$e_{\theta}^*$', 'Interpreter', 'latex')
    legend boxoff
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',18)
    axis([0 Tmax 0 0.6])
    box
    subplot(3,1,2)
    hold on
    e_L_line=plot(TSample,e_L(:,1), 'r','LineWidth',3);
    e_L_star_line=yline(compactness_thresh,'r--','LineWidth',2);
    axis([0 Tmax 0 0.9])
    yticks([0 .3 .6 .9])
    legend('$\bar{e}_{L}$', '$e_{L}^*$', 'Interpreter', 'latex')
    legend boxoff
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',18)
    box
    subplot(3,1,3)
    plot([0 Tmax/3 Tmax/3 2*Tmax/3 2*Tmax/3 Tmax], [4 4 6 6 4 4], 'k', 'LineWidth', 2)
    axis([0 Tmax 3 7])
    yticks([4 6])
    box
    legend('$L$', 'Interpreter', 'latex')
    legend boxoff
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',18)
    set(gcf,'Position',[100 100 500 500])
    xlabel('t')
end


end

