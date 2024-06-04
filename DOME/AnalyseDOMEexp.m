clear
close all

%% Set the folder that you want to analyse 


% Add the path to your data
data_folder = 'C:\Users\david\OneDrive - Università di Napoli Federico II\Research\Data\DOME\Euglena_switch_10\combo5'; % 10 s Full experiment  
% data_folder = ['C:\Users\david\OneDrive - Università di Napoli Federico II\Research\Data\DOME\Euglena_255_ON\combo'];
% data_folder = ['C:\Users\david\OneDrive - Università di Napoli Federico II\Research\Data\DOME\Euglena_OFF\combo'];


% FOR FAST SWITCH BETWEEN FOLDERS
% data_folder = '/Volumes/DOMEPEN/Experiments/comparisons/Euglena_OFF/combo'; % off combo
% data_folder = '/Volumes/DOMEPEN/Experiments/comparisons/Euglena_switch_5/combo'; % switch5s combo
% data_folder = '/Volumes/DOMEPEN/Experiments/comparisons/Euglena_255_ON/combo'; % on255 combo
% data_folder = '/Volumes/DOMEPEN/Experiments/2023_06_15_Euglena_7/tracking_2023_10_16'; % switch10s
% data_folder = '/Volumes/DOMEPEN/Experiments/comparisons/Euglena_switch_10/combo5'; % switch10s combo
% data_folder = '/Volumes/DOMEPEN/Experiments/2023_06_26_Euglena_19/tracking_2023_10_16'; % on255



%% Load data and set experiment parameters


deltaT = 0.5;
dT = 0.01;

current_folder = fileparts(which('Launcher'));
addpath(genpath(current_folder));


% Load longitudinal, angular velocities and the light inputs in time
speed  = load(fullfile(data_folder,'speeds_smooth.txt'));
omega  = load(fullfile(data_folder,'ang_vel_smooth.txt'));
inputs = load(fullfile(data_folder,'inputs.txt'));


plot_data = true;                                   %Flag to plot the experimental data
stat_an = false;
N = size(speed,2);                                  %Number of Agents
timeInstants = [0:size(speed,1)-1] * deltaT;        %Time vecor
agents = [0:N-1]';                                  %Agents ID
u=inputs(:,1)/255;                                  %Normalized input (0,1)
u_dot_BE = [0;diff(u)]/deltaT;                      %Derivative of the control input
% u_dot_grad = gradient(u)/deltaT;
u_dot = u_dot_BE;                                   %Standardise notation
u_dotn = min(u_dot,0);                              %Study step-down response
u_dotp = max(u_dot,0);                              %Study step-up response
u_matrix = [u, u_dot];                              %Input vecotrization



%% Statistical analysis of the signal (All points)

if stat_an


    switchn_t = find(u_dotn);                                %Find all the switch on times
    t_a = 5;                                                 %Time after the switch I want to analyse
    s_aftern = [];                                           %Speed after the switch
    s_beforen =[];                                           %Speed before the switch
    w_aftern = [];                                           %Angular velocity after the switch
    w_beforen =[];                                           %Angular velocity before the switch

    switchp_t = find(u_dotp);                                %Find all the switch on times
    s_afterp = [];                                           %Speed after the switch
    s_beforep =[];                                           %Speed before the switch
    w_afterp = [];                                           %Angular velocity after the switch
    w_beforep =[];                                           %Angular velocity before the switch

    for i=1:length(switchn_t)-1
        s_aftern  = [s_aftern;  speed(switchn_t(i):switchn_t(i)+round(t_a/deltaT),:)];
        s_beforen = [s_beforen; speed(switchn_t(i)-round(t_a/deltaT):switchn_t(i),:)];
        w_aftern  = [w_aftern;  abs(omega(switchn_t(i):switchn_t(i)+round(t_a/deltaT),:))];
        w_beforen = [w_beforen; abs(omega(switchn_t(i)-round(t_a/deltaT):switchn_t(i),:))];
    end


    for i=1:length(switchp_t)-1
        s_afterp  = [s_afterp;  speed(switchp_t(i):switchp_t(i)+round(t_a/deltaT),:)];
        s_beforep = [s_beforep; speed(switchp_t(i)-round(t_a/deltaT):switchp_t(i),:)];
        w_afterp  = [w_afterp;  abs(omega(switchp_t(i):switchp_t(i)+round(t_a/deltaT),:))];
        w_beforep = [w_beforep; abs(omega(switchp_t(i)-round(t_a/deltaT):switchp_t(i),:))];
    end

    s_beforen_avg = mean(s_beforen,1,"omitnan");
    s_aftern_avg = mean(s_aftern,1,"omitnan");

    w_beforen_avg = mean(w_beforen,1,"omitnan");
    w_aftern_avg = mean(w_aftern,1,"omitnan");


    s_beforep_avg = mean(s_beforep,1,"omitnan");
    s_afterp_avg = mean(s_afterp,1,"omitnan");

    w_beforep_avg = mean(w_beforep,1,"omitnan");
    w_afterp_avg = mean(w_afterp,1,"omitnan");


    %Parameters of the plots
    k_est = true;                                   %Show continous densities
    n_bins = 16;                                    %Bins of the Histogram
    Pix_SS = get(0,'screensize');                   %Get screen n of pixels
    sc_w = Pix_SS(3);                               %Get screen width
    sc_h = Pix_SS(4);                               %Get screen heigt
    offs = 50;                                      %Offset to compensate windows bar
    sc_h = sc_h-offs;
    avg_an = true;


    %Analysis of the average in time
    if avg_an

        s_beforen = s_beforen_avg;
        s_aftern = s_aftern_avg;

        w_beforen = w_beforen_avg;
        w_aftern = w_aftern_avg;


        s_beforep = s_beforep_avg;
        s_afterp = s_afterp_avg;

        w_beforep = w_beforep_avg;
        w_afterp = w_afterp_avg;

    end






    %%%%%%%%%%%%%%%%%%%%%%%%%%%% TESTING ALPHA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Comparing the distributions of the data n seconds before the switch ON and
    %n seconds before the switch OFF

    %FIRST pair : v, Second pair: w


    figure("Name","Speed alpha","Position",[0 offs sc_w/4 sc_h/3]);


    edges = linspace(min([s_beforep(:);s_beforen(:)]),max([s_beforep(:);s_beforen(:)]),n_bins+1);
    myHistogram(s_beforep(:),edges,k_est);
    hold on;
    myHistogram(s_beforen(:),edges,k_est);
    xlabel('Speed(px/s)');
    % ylabel('Probability');
    % legend(sprintf("OFF (%ds before switch)",t_a),sprintf("ON (%ds before switch)",t_a))
    set(gca,'FontSize',14);


    figure("Name","Speed alpha (BoxPlot)","Position",[sc_w/4 offs sc_w/4 sc_h/3]);

    myboxplot({s_beforep(:),s_beforen(:)},true);
    xticklabels({sprintf("OFF (%ds)",t_a),sprintf("ON (%ds)",t_a)})
    ylabel('Speed(px/s)');
    set(gca,'FontSize',14);



    figure("Name","Omega alpha (Histogram)","Position",[2*sc_w/4 offs sc_w/4 sc_h/3]);

    edges = linspace(min([w_beforep(:);w_beforen(:)]),max([w_beforep(:);w_beforen(:)]),n_bins+1);
    myHistogram(w_beforep(:),edges,k_est);
    hold on;
    myHistogram(w_beforen(:),edges,k_est);
    xlabel('\omega(rad/s)');
    % ylabel('Probability');
    % legend(sprintf("OFF (%ds before switch)",t_a),sprintf("ON (%ds before switch)",t_a))
    set(gca,'FontSize',14);


    figure("Name","Omega alpha (BoxPlot)","Position",[3*sc_w/4 offs sc_w/4 sc_h/3]);

    myboxplot({w_beforep(:),w_beforen(:)},true);
    xticklabels({sprintf("OFF (%ds)",t_a),sprintf("ON (%ds)",t_a)})
    ylabel('\omega(rad/s)');
    set(gca,'FontSize',14);




    %%%%%%%%%%%%%%%%%%%%%%%%%%%% TESTING BETA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Comparing the distributions of the data n seconds before the switch ON and
    %n seconds before the switch OFF

    %FIRST pair : v, Second pair: w


    figure("Name","Speed beta (Histogram)","Position",[0 offs+sc_h/3 sc_w/4 sc_h/3]);


    edges = linspace(min([s_afterp(:);s_beforen(:)]),max([s_afterp(:);s_beforen(:)]),n_bins+1);
    myHistogram(s_afterp(:),edges,k_est);
    hold on;
    myHistogram(s_beforen(:),edges,k_est);
    xlabel('Speed(px/s)');
    % ylabel('Probability');
    % legend(sprintf("ON (%ds after switch)",t_a),sprintf("ON (%ds before switch)",t_a))
    set(gca,'FontSize',14);


    figure("Name","Speed beta (BoxPlot)","Position",[sc_w/4 offs+sc_h/3 sc_w/4 sc_h/3]);

    myboxplot({s_afterp(:),s_beforen(:)},true);
    xticklabels({sprintf("ON as (%ds)",t_a),sprintf("ON bs (%ds)",t_a)})
    ylabel('Speed(px/s)');
    set(gca,'FontSize',14);



    figure("Name","Omega beta (Histogram)","Position",[2*sc_w/4 offs+sc_h/3 sc_w/4 sc_h/3]);

    edges = linspace(min([w_afterp(:);w_beforen(:)]),max([w_afterp(:);w_beforen(:)]),n_bins+1);
    myHistogram(w_afterp(:),edges,k_est);
    hold on;
    myHistogram(w_beforen(:),edges,k_est);
    xlabel('\omega(rad/s)');
    % ylabel('Probability');
    % legend(sprintf("ON (%ds after switch)",t_a),sprintf("ON (%ds before switch)",t_a))
    set(gca,'FontSize',14);


    figure("Name","Omega alpha (BoxPlot)","Position",[3*sc_w/4 offs+sc_h/3 sc_w/4 sc_h/3]);

    myboxplot({w_afterp(:),w_beforen(:)},true);
    xticklabels({sprintf("ON as (%ds)",t_a),sprintf("ON bs (%ds)",t_a)})
    ylabel('\omega(rad/s)');
    set(gca,'FontSize',14);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%% TESTING GAMMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Comparing the distributions of the data n seconds before the switch ON and
    %n seconds after the switch OFF

    %FIRST pair : v, Second pair: w


    figure("Name","Speed \gamma (Histogram)","Position",[0 offs+2*sc_h/3 sc_w/4 sc_h/3]);


    edges = linspace(min([s_beforep(:);s_aftern(:)]),max([s_beforep(:);s_aftern(:)]),n_bins+1);
    myHistogram(s_aftern(:),edges,k_est);
    hold on;
    myHistogram(s_beforep(:),edges,k_est);
    xlabel('Speed(px/s)');
    % ylabel('Probability');
    % legend(sprintf("OFF (%ds after switch)",t_a),sprintf("OFF (%ds before switch)",t_a))
    set(gca,'FontSize',14);


    figure("Name","Speed gamma (BoxPlot)","Position",[sc_w/4 offs+2*sc_h/3 sc_w/4 sc_h/3]);

    myboxplot({s_aftern(:),s_beforep(:)},true);
    xticklabels({sprintf("OFF an (%ds)",t_a),sprintf("OFF bp (%ds)",t_a)})
    ylabel('Speed(px/s)');
    set(gca,'FontSize',14);



    figure("Name","\omega \beta (Histogram)","Position",[2*sc_w/4 offs+2*sc_h/3 sc_w/4 sc_h/3]);

    edges = linspace(min([w_beforen(:);w_afterp(:)]),max([w_beforen(:);w_afterp(:)]),n_bins+1);
    myHistogram(w_aftern(:),edges,k_est);
    hold on;
    myHistogram(w_beforep(:),edges,k_est);
    xlabel('\omega(rad/s)');
    % ylabel('Probability');
    % legend(sprintf("OFF (%ds after switch)",t_a),sprintf("OFF (%ds before switch)",t_a))
    set(gca,'FontSize',14);


    figure("Name","Omega alpha (BoxPlot)","Position",[3*sc_w/4 offs+2*sc_h/3 sc_w/4 sc_h/3]);

    myboxplot({w_aftern(:),w_beforep(:)},true);
    xticklabels({sprintf("OFF ap (%ds)",t_a),sprintf("ON bn (%ds)",t_a)})
    ylabel('\omega(rad/s)');
    set(gca,'FontSize',14);


    pause;
    close all;

end


%% PLOTS
%% 

if plot_data


    outputDir = fullfile(data_folder,'plots');
    if ~exist(outputDir,'dir'); mkdir(outputDir); end

    figure % SCATTER PLOT - SPEED and ANGULAR VELOCITY - MEAN OVER TIME
    x_to_plot=[];
    y_to_plot=[];
    grouping=[];
    s=scatterhist(mean(speed,1,'omitnan'), mean(abs(omega),1,'omitnan'), 'Location','NorthEast','Direction','out' ,'Color','bk','Kernel','on');
    xlabel(s,'mean $v$ [$\mu$m/s]','Interpreter','Latex','FontSize',16)
    ylabel(s,'mean $|\omega|$ [rad/s]','Interpreter','Latex','FontSize',16)
    s(1).YAxisLocation = 'left';
    s(1).XAxisLocation = 'bottom';
    %set(get(gca,'children'),'filled',true)
    s(2).Position = [0.1    0.82   0.7    0.125];
    s(3).Position = [0.82   0.1    0.125    0.7];
    s(1).Position(3) = 0.7;
    s(1).Position(4) = 0.7;
    ylim([0,3])
    xlim([0,120])
    if outputDir
        saveas(gcf,fullfile(outputDir, 'scatter_meanOnTime'))
        saveas(gcf,fullfile(outputDir, 'scatter_meanOnTime'),'png')
    end

    figure % TIME PLOT - SPEED and ANGULAR VELOCITY
    subplot(2,1,1)
    xlim([0,max(timeInstants)])
    ylim([0,120])
    if isvarname('u')
        highlightInputs(timeInstants, u, 'r', 0.25)
    end
    l1=plotWithShade(timeInstants, median(speed,2,'omitnan'), quantile(speed, 0.1, 2), quantile(speed, 0.9, 2), 'b', 0.3);
    xlabel('$t$ [s]','Interpreter','Latex','FontSize',16)
    ylabel('$v$ [$\mu$m/s]','Interpreter','Latex','FontSize',16)
    rng=ylim;
    box on
    subplot(2,1,2)
    xlim([0,max(timeInstants)])
    ylim([0,2])
    if isvarname('u')
        highlightInputs(timeInstants, u, 'r', 0.25)
    end
    l1=plotWithShade(timeInstants, median(abs(omega),2,'omitnan'), quantile(abs(omega), 0.1, 2), quantile(abs(omega), 0.9, 2), 'b', 0.3);
    xlabel('$t$ [s]','Interpreter','Latex','FontSize',16)
    ylabel('$|\omega|$ [rad/s]','Interpreter','Latex','FontSize',16)
    rng=ylim;
    box on
    if outputDir
        saveas(gcf,fullfile(outputDir, 'time_plot'))
        saveas(gcf,fullfile(outputDir, 'time_plot'),'png')
    end

    figure % BOX PLOT - SPEED and ANGULAR VELOCITY - Mean over time
    hold on
    set(gca,'FontSize',12)
    yline(0)
    myboxplot({mean(omega,1,'omitnan')}, true, 3, {'b'})%, [0,0.4470,0.7410])
    % myboxplot({omega(:)}, true, 3, {'b'})%, [0,0.4470,0.7410])
    xticks([])
    ylabel('$\omega$ [rad/s]','Interpreter','Latex','FontSize',16)
    set(gcf,'position',[300,300,300,420])
    if outputDir
        saveas(gcf,fullfile(outputDir, 'ang_vel_boxplot_meanovertime'))
        saveas(gcf,fullfile(outputDir, 'ang_vel_boxplot_meanovertime'),'png')
    end

    
    figure % BOX PLOT - SPEED and ANGULAR VELOCITY - Mean over time
    hold on
    set(gca,'FontSize',12)
    yline(0)
    myboxplot({mean(speed,1,'omitnan')}, true, 3, {'b'})%, [0,0.4470,0.7410])
    % myboxplot({omega(:)}, true, 3, {'b'})%, [0,0.4470,0.7410])
    xticks([])
    ylabel('$v$ [px/s]','Interpreter','Latex','FontSize',16)
    set(gcf,'position',[300,300,300,420])


    figure
    hold on
    yline(0)
    l1=plotWithShade(timeInstants, median(omega,2,'omitnan'), quantile(omega, 0.1, 2), quantile(omega, 0.9, 2), 'b', 0.3);
    yline(median(median(omega,2,'omitnan')),'b')
    if isvarname('u')
        highlightInputs(timeInstants, u, 'r', 0.25)
    end


    % figure % inputs
    % subplot(3,1,1)
    % hold on
    % plot(timeInstants,u,'--')
    % plot(t_sim,u_sim(:,1),'--')
    % legend({'experiment','simulation'})
    % subplot(3,1,2)
    % hold on
    % plot(timeInstants,u_dot,'--')
    % plot(t_sim,u_sim(:,2),'--')
    % legend({'experiment','simulation'})
    % subplot(3,1,3)
    % hold on
    % plot(timeInstants,cumtrapz(u_dot)*deltaT,'--')
    % plot(t_sim,cumtrapz(u_sim(:,2))*dT,'--')
    % legend({'experiment','simulation'})


end