clear
close all

%% Set the folder that you want to analyse


% Add the path to your data
% experiments_folder="C:\Users\david\OneDrive - Università di Napoli Federico II\Research\Data\DOME\";    % DAVIDE
experiments_folder="/Volumes/DOMEPEN/Experiments/comparisons";                                          % ANDREA

% experiments_names=[fullfile("Euglena_75_ON","combo"),fullfile("Euglena_150_ON","combo"),fullfile("Euglena_255_ON","combo")];
% experiments_names=[fullfile("Euglena_255_ON","combo")];
% experiments_names=[fullfile("Euglena_OFF","combo")];

% experiments_names=[fullfile("Volvox_75_ON","combo"),fullfile("Volvox_150_ON","combo"),fullfile("Volvox_255_ON","combo")];
experiments_names=[fullfile("Volvox_255_ON","combo")];
experiments_names=[fullfile("Volvox_switch_10","combo5")];
% experiments_names=[fullfile("Volvox_OFF","combo")];


plot_data = true;                                   % Plot the experimental data
stat_an = false;                                      % Plot statistical analysis of light response
avg_an = true;                                       % Use average values (over times) for stat analysis

t_a = 12;                                            % Time window after the switch I want to analyse
t_b = 12;                                            % Time window before the switch I want to analyse


deltaT = 0.5;

speed_lim = 200;
omega_lim = 1.75;


% allocate vars
s_afterp = cell(1,length(experiments_names));                                           %Speed after the switch
s_beforep =cell(1,length(experiments_names));                                           %Speed before the switch
w_afterp = cell(1,length(experiments_names));                                           %Angular velocity after the switch
w_beforep =cell(1,length(experiments_names));                                           %Angular velocity before the switch
s_aftern = cell(1,length(experiments_names));                                           %Speed after the switch
s_beforen =cell(1,length(experiments_names));                                           %Speed before the switch
w_aftern = cell(1,length(experiments_names));                                           %Angular velocity after the switch
w_beforen =cell(1,length(experiments_names));                                           %Angular velocity before the switch



for exp=1:length(experiments_names)
    
    %% Load data and set experiment parameters
    
    
    % current_folder = fileparts(which('Launcher'));
    % addpath(genpath(current_folder));
    
    data_folder = fullfile(experiments_folder,experiments_names(exp));
    % Load longitudinal, angular velocities and the light inputs in time
    speed  = load(fullfile(data_folder,'speeds_smooth.txt'));
    omega  = load(fullfile(data_folder,'ang_vel_smooth.txt'));
    inputs = load(fullfile(data_folder,'inputs.txt'));
    
    
    % speed = speed./median(median(speed(1:60,:),2,'omitnan'),'omitnan');
    % omega = abs(omega)./median(median(abs(omega(1:60,:)),2,'omitnan'),'omitnan');
    
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
        switchp_t = find(u_dotp);                                %Find all the switch on times
        
        for i=1:length(switchn_t)
            s_aftern{exp}  = [s_aftern{exp};  speed(switchn_t(i):switchn_t(i)+round(t_a/deltaT),:)];
            s_beforen{exp} = [s_beforen{exp}; speed(switchn_t(i)-round(t_b/deltaT):switchn_t(i),:)];
            w_aftern{exp}  = [w_aftern{exp};  abs(omega(switchn_t(i):switchn_t(i)+round(t_a/deltaT),:))];
            w_beforen{exp} = [w_beforen{exp}; abs(omega(switchn_t(i)-round(t_b/deltaT):switchn_t(i),:))];
        end
        
        
        for i=1:length(switchp_t)
            s_afterp{exp}  = [s_afterp{exp};  speed(switchp_t(i):switchp_t(i)+round(t_a/deltaT),:)];
            s_beforep{exp} = [s_beforep{exp}; speed(switchp_t(i)-round(t_b/deltaT):switchp_t(i),:)];
            w_afterp{exp}  = [w_afterp{exp};  abs(omega(switchp_t(i):switchp_t(i)+round(t_a/deltaT),:))];
            w_beforep{exp} = [w_beforep{exp}; abs(omega(switchp_t(i)-round(t_b/deltaT):switchp_t(i),:))];
        end
        
        s_beforen_avg = mean(s_beforen{exp},1,"omitnan");
        s_aftern_avg = mean(s_aftern{exp},1,"omitnan");
        
        w_beforen_avg = mean(w_beforen{exp},1,"omitnan");
        w_aftern_avg = mean(w_aftern{exp},1,"omitnan");
        
        
        s_beforep_avg = mean(s_beforep{exp},1,"omitnan");
        s_afterp_avg = mean(s_afterp{exp},1,"omitnan");
        
        w_beforep_avg = mean(w_beforep{exp},1,"omitnan");
        w_afterp_avg = mean(w_afterp{exp},1,"omitnan");
        
        
        %Parameters of the plots
        k_est = true;                                   %Show continous densities
        n_bins = 16;                                    %Bins of the Histogram
        Pix_SS = get(0,'screensize');                   %Get screen n of pixels
        sc_w = Pix_SS(3);                               %Get screen width
        sc_h = Pix_SS(4);                               %Get screen heigt
        offs = 50;                                      %Offset to compensate windows bar
        sc_h = sc_h-offs;
        
        
        %Analysis of the average in time
        if avg_an
            
            s_beforen{exp} = s_beforen_avg;
            s_aftern{exp} = s_aftern_avg;
            
            w_beforen{exp} = w_beforen_avg;
            w_aftern{exp} = w_aftern_avg;
            
            
            s_beforep{exp} = s_beforep_avg;
            s_afterp{exp} = s_afterp_avg;
            
            w_beforep{exp} = w_beforep_avg;
            w_afterp{exp} = w_afterp_avg;
            
        end
        
        
        
        
        
        outputDir = fullfile(data_folder,'plots');
        if ~exist(outputDir,'dir'); mkdir(outputDir); end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% TESTING ALPHA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Comparing the distributions of the data n seconds before the switch ON and
        %n seconds before the switch OFF
        
        %FIRST pair : v, Second pair: w
        
        
        figure("Name","Speed alpha (Histogram)","Position",[0 offs sc_w/4 sc_h/3]);
        edges = linspace(min([s_beforep{exp},s_beforen{exp}],[],"omitnan"),max([s_beforep{exp},s_beforen{exp}],[],"omitnan"),n_bins+1);
        myHistogram(s_beforep{exp},edges,k_est);
        hold on;
        myHistogram(s_beforen{exp},edges,k_est);
        xlabel('Speed(px/s)');
        xlim([0,speed_lim]);
        % ylabel('Probability');
        % legend(sprintf("OFF (%ds before switch)",t_b),sprintf("ON (%ds before switch)",t_b))
        set(gca,'FontSize',14);
        
        
        figure("Name","Speed alpha (BoxPlot)","Position",[sc_w/4 offs sc_w/4 sc_h/3]);
        myboxplot({s_beforep{exp},s_beforen{exp}},true);
        xticklabels({sprintf("OFF (%ds)",t_b),sprintf("ON (%ds)",t_b)})
        ylabel('Speed(px/s)');
        ylim([0,speed_lim]);
        set(gca,'FontSize',14);
        if isfolder(outputDir)
            saveas(gcf,fullfile(outputDir, 'boxplot_alpha_speed'))
            saveas(gcf,fullfile(outputDir, 'boxplot_alpha_speed'),'png')
        end
        
        
        figure("Name","Omega alpha (Histogram)","Position",[2*sc_w/4 offs sc_w/4 sc_h/3]);
        edges = linspace(min([w_beforep{exp},w_beforen{exp}],[],"omitnan"),max([w_beforep{exp},w_beforen{exp}],[],"omitnan"),n_bins+1);
        myHistogram(w_beforep{exp},edges,k_est);
        hold on;
        myHistogram(w_beforen{exp},edges,k_est);
        xlabel('\omega(rad/s)');
        xlim([0,omega_lim]);
        % ylabel('Probability');
        % legend(sprintf("OFF (%ds before switch)",t_a),sprintf("ON (%ds before switch)",t_a))
        set(gca,'FontSize',14);
        
        
        figure("Name","Omega alpha (BoxPlot)","Position",[3*sc_w/4 offs sc_w/4 sc_h/3]);
        myboxplot({w_beforep{exp},w_beforen{exp}},true);
        xticklabels({sprintf("OFF (%ds)",t_b),sprintf("ON (%ds)",t_b)})
        ylabel('\omega(rad/s)');
        ylim([0,omega_lim]);
        set(gca,'FontSize',14);
        if isfolder(outputDir)
            saveas(gcf,fullfile(outputDir, 'boxplot_alpha_omega'))
            saveas(gcf,fullfile(outputDir, 'boxplot_alpha_omega'),'png')
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% TESTING BETA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Comparing the distributions of the data n seconds before the switch ON and
        %n seconds before the switch OFF
        
        %FIRST pair : v, Second pair: w
        
        
        figure("Name","Speed beta (Histogram)","Position",[0 offs+sc_h/3 sc_w/4 sc_h/3]);
        edges = linspace(min([s_afterp{exp},s_beforen{exp}],[],"omitnan"),max([s_afterp{exp},s_beforen{exp}],[],"omitnan"),n_bins+1);
        myHistogram(s_afterp{exp},edges,k_est);
        hold on;
        myHistogram(s_beforen{exp},edges,k_est);
        xlabel('Speed(px/s)');
        xlim([0,speed_lim]);
        % ylabel('Probability');
        % legend(sprintf("ON (%ds after switch)",t_a),sprintf("ON (%ds before switch)",t_a))
        set(gca,'FontSize',14);
        
        
        figure("Name","Speed beta (BoxPlot)","Position",[sc_w/4 offs+sc_h/3 sc_w/4 sc_h/3]);
        myboxplot({s_afterp{exp},s_beforen{exp}},true);
        xticklabels({sprintf("ON as (%ds)",t_a),sprintf("ON bs (%ds)",t_b)})
        ylabel('Speed(px/s)');
        ylim([0,speed_lim]);
        set(gca,'FontSize',14);
        if isfolder(outputDir)
            saveas(gcf,fullfile(outputDir, 'boxplot_beta_speed'))
            saveas(gcf,fullfile(outputDir, 'boxplot_beta_speed'),'png')
        end
        
        
        figure("Name","Omega beta (Histogram)","Position",[2*sc_w/4 offs+sc_h/3 sc_w/4 sc_h/3]);
        edges = linspace(min([w_afterp{exp},w_beforen{exp}],[],"omitnan"),max([w_afterp{exp},w_beforen{exp}],[],"omitnan"),n_bins+1);
        myHistogram(w_afterp{exp},edges,k_est);
        hold on;
        myHistogram(w_beforen{exp},edges,k_est);
        xlabel('\omega(rad/s)');
        xlim([0,omega_lim]);
        % ylabel('Probability');
        % legend(sprintf("ON (%ds after switch)",t_a),sprintf("ON (%ds before switch)",t_a))
        set(gca,'FontSize',14);
        
        
        figure("Name","Omega beta (BoxPlot)","Position",[3*sc_w/4 offs+sc_h/3 sc_w/4 sc_h/3]);
        myboxplot({w_afterp{exp},w_beforen{exp}},true);
        xticklabels({sprintf("ON as (%ds)",t_a),sprintf("ON bs (%ds)",t_b)})
        ylabel('\omega(rad/s)');
        ylim([0,omega_lim]);
        set(gca,'FontSize',14);
        if isfolder(outputDir)
            saveas(gcf,fullfile(outputDir, 'boxplot_beta_omega'))
            saveas(gcf,fullfile(outputDir, 'boxplot_beta_omega'),'png')
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% TESTING GAMMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Comparing the distributions of the data n seconds before the switch ON and
        %n seconds after the switch OFF
        
        %FIRST pair : v, Second pair: w
        
        
        figure("Name","Speed gamma (Histogram)","Position",[0 offs+2*sc_h/3 sc_w/4 sc_h/3]);
        edges = linspace(min([s_beforep{exp},s_aftern{exp}],[],"omitnan"),max([s_beforep{exp},s_aftern{exp}],[],"omitnan"),n_bins+1);
        myHistogram(s_aftern{exp},edges,k_est);
        hold on;
        myHistogram(s_beforep{exp},edges,k_est);
        xlabel('Speed(px/s)');
        xlim([0,speed_lim]);
        % ylabel('Probability');
        % legend(sprintf("OFF (%ds after switch)",t_a),sprintf("OFF (%ds before switch)",t_a))
        set(gca,'FontSize',14);
        
        
        figure("Name","Speed gamma (BoxPlot)","Position",[sc_w/4 offs+2*sc_h/3 sc_w/4 sc_h/3]);
        myboxplot({s_aftern{exp},s_beforep{exp}},true);
        xticklabels({sprintf("OFF an (%ds)",t_a),sprintf("OFF bp (%ds)",t_b)})
        ylabel('Speed(px/s)');
        ylim([0,speed_lim]);
        set(gca,'FontSize',14);
        if isfolder(outputDir)
            saveas(gcf,fullfile(outputDir, 'boxplot_gamma_speed'))
            saveas(gcf,fullfile(outputDir, 'boxplot_gamma_speed'),'png')
        end
        
        
        figure("Name","Omega gamma (Histogram)","Position",[2*sc_w/4 offs+2*sc_h/3 sc_w/4 sc_h/3]);
        edges = linspace(min([w_beforen{exp},w_afterp{exp}],[],"omitnan"),max([w_beforen{exp},w_afterp{exp}],[],"omitnan"),n_bins+1);
        myHistogram(w_aftern{exp},edges,k_est);
        hold on;
        myHistogram(w_beforep{exp},edges,k_est);
        xlabel('\omega(rad/s)');
        xlim([0,omega_lim]);
        % ylabel('Probability');
        % legend(sprintf("OFF (%ds after switch)",t_a),sprintf("OFF (%ds before switch)",t_a))
        set(gca,'FontSize',14);
        
        
        figure("Name","Omega gamma (BoxPlot)","Position",[3*sc_w/4 offs+2*sc_h/3 sc_w/4 sc_h/3]);
        myboxplot({w_aftern{exp},w_beforep{exp}},true);
        xticklabels({sprintf("OFF ap (%ds)",t_a),sprintf("ON bn (%ds)",t_b)})
        ylabel('\omega(rad/s)');
        ylim([0,omega_lim]);
        set(gca,'FontSize',14);
        if isfolder(outputDir)
            saveas(gcf,fullfile(outputDir, 'boxplot_gamma_omega'))
            saveas(gcf,fullfile(outputDir, 'boxplot_gamma_omega'),'png')
        end
        
        
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
        ylim([0,omega_lim])
        xlim([0,speed_lim])
        if isfolder(outputDir)
            saveas(gcf,fullfile(outputDir, 'scatter_meanOnTime'))
            saveas(gcf,fullfile(outputDir, 'scatter_meanOnTime'),'png')
        end
        
        figure % TIME PLOT - SPEED and ANGULAR VELOCITY
        subplot(2,1,1)
        hold on
        xlim([0,max(timeInstants)])
        ylim([0,speed_lim])
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
        ylim([0,omega_lim])
        if isvarname('u')
            highlightInputs(timeInstants, u, 'r', 0.25)
        end
        l1=plotWithShade(timeInstants, median(abs(omega),2,'omitnan'), quantile(abs(omega), 0.1, 2), quantile(abs(omega), 0.9, 2), 'b', 0.3);
        xlabel('$t$ [s]','Interpreter','Latex','FontSize',16)
        ylabel('$|\omega|$ [rad/s]','Interpreter','Latex','FontSize',16)
        rng=ylim;
        box on
        if isfolder(outputDir)
            saveas(gcf,fullfile(outputDir, 'time_plot'))
            saveas(gcf,fullfile(outputDir, 'time_plot'),'png')
        end
        
        figure % BOX PLOT - ANGULAR VELOCITY - Mean over time
        hold on
        set(gca,'FontSize',12)
        yline(0)
        myboxplot({mean(omega,1,'omitnan')}, true, 3, {'b'})%, [0,0.4470,0.7410])
        % myboxplot({omega(:)}, true, 3, {'b'})%, [0,0.4470,0.7410])
        xticks([])
        ylim([-omega_lim,omega_lim])
        ylabel('$\omega$ [rad/s]','Interpreter','Latex','FontSize',16)
        set(gcf,'position',[300,300,300,420])
        if isfolder(outputDir)
            saveas(gcf,fullfile(outputDir, 'ang_vel_boxplot_meanovertime'))
            saveas(gcf,fullfile(outputDir, 'ang_vel_boxplot_meanovertime'),'png')
        end
        
        
        figure % BOX PLOT - SPEED - Mean over time
        hold on
        set(gca,'FontSize',12)
        yline(0)
        myboxplot({mean(speed,1,'omitnan')}, true, 3, {'b'})%, [0,0.4470,0.7410])
        % myboxplot({omega(:)}, true, 3, {'b'})%, [0,0.4470,0.7410])
        xticks([])
        ylabel('$v$ [px/s]','Interpreter','Latex','FontSize',16)
        set(gcf,'position',[300,300,300,420])
        if isfolder(outputDir)
            saveas(gcf,fullfile(outputDir, 'speed_boxplot_meanovertime'))
            saveas(gcf,fullfile(outputDir, 'speed_boxplot_meanovertime'),'png')
        end
        
        %         figure
        %         hold on
        %         yline(0)
        %         l1=plotWithShade(timeInstants, median(abs(omega),2,'omitnan'), quantile(abs(omega), 0.1, 2), quantile(abs(omega), 0.9, 2), 'b', 0.3);
        %         yline(median(median(omega,2,'omitnan')),'b')
        %         if isvarname('u')
        %             highlightInputs(timeInstants, u, 'r', 0.25)
        %         end
        
        %         pause;
        %         close all;
    end
    
    
end



% compare across experiments
if length(experiments_names)>1
    figure();
    for i=1:length(experiments_names)
        s_afterp{i}= s_afterp{i}/median(s_beforen{i},"omitnan");
        sm(i) = median(s_afterp{i},"omitnan");
        bp(i) = quantile(s_afterp{i},0.25);
        up(i) = quantile(s_afterp{i},0.75);
    end
    plotWithShade([0 75/255 150/255 255/255],[1 sm],[1 bp],[1 up],[1 0 0],0.3);
    title('BetaS')
    
    figure();
    for i=1:length(experiments_names)
        s_beforen{i}= s_beforen{i}/median(s_beforep{i},"omitnan");
        sm(i) = median(s_beforen{i},"omitnan");
        bp(i) = quantile(s_beforen{i},0.25);
        up(i) = quantile(s_beforen{i},0.75);
    end
    plotWithShade([0 75/255 150/255 255/255],[1 sm],[1 bp],[1 up],[1 0 0],0.3);
    title('AlphaS')
    
    figure();
    for i=1:length(experiments_names)
        w_afterp{i}= w_afterp{i}/median(w_beforen{i},"omitnan");
        sm(i) = median(w_afterp{i},"omitnan");
        bp(i) = quantile(w_afterp{i},0.25);
        up(i) = quantile(w_afterp{i},0.75);
    end
    
    plotWithShade([0 75/255 150/255 255/255],[1 sm],[1 bp],[1 up],[1 0 0],0.3);
    title('BetaW')
    
    figure();
    for i=1:length(experiments_names)
        w_beforen{i}= w_beforen{i}/median(w_beforep{i},"omitnan");
        sm(i) = median(w_beforen{i},"omitnan");
        bp(i) = quantile(w_beforen{i},0.25);
        up(i) = quantile(w_beforen{i},0.75);
        
    end
    
    plotWithShade([0 75/255 150/255 255/255],[1 sm],[1 bp],[1 up],[1 0 0],0.3);
    title('AlphaW')
end
