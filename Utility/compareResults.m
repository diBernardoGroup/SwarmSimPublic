%% COMPARE EXPERIMENTS

function [metrics_of_interest] = compareResults(experiment_paths, outputDir, make_plots, Render)

arguments
    experiment_paths
    outputDir = ''
    make_plots = true
    Render = struct()
end

%% Add subfolders to the Matlab path
current_folder = fileparts(which('defaultParam'));
addpath(genpath(current_folder));

number_of_exp = length(experiment_paths);
experiments={};

%% Load data
timeInstants = 0:0.5:180;
for i=1:number_of_exp
    exp = experiment_paths{i};
    data_path = fullfile(exp,'data.mat');
    if isfile(data_path) % load data from simulation
        d = load(data_path);
        omega = d.omega';
        speed = d.speed';
        if isfield(d,'timeInstants')
            timeInstants = d.timeInstants;
        else
            timeInstants = d.Simulation.timeInstants;
        end
        if i==1
            u = d.u;
        else
            assert(all(u == d.u),'Experiments must have the same inputs')
        end
    else                % load data from experiment
        if ~isfile(fullfile(exp,'speeds_smooth.txt'))
            tracking = getLastTracking(exp);
            exp = fullfile(exp,tracking);
        end
        speed = load(fullfile(exp,'speeds_smooth.txt'))';
        omega = load(fullfile(exp,'ang_vel_smooth.txt'))';
        inputs=load(fullfile(exp,'inputs.txt'));
        if i==1
            u=inputs(:,1)/255;              %select blue channel and scale in [0,1]
        else
            assert(all(u == inputs(:,1)/255),'Experiments must have the same inputs')
        end
    end
    
    data = table(speed, omega);
    data(all(isnan(data.speed),2),:) = [];
    experiments{i}=data;
end

if number_of_exp>2
    nmse_speed = nan(number_of_exp);
    nmse_omega = nan(number_of_exp);
    for i=1:number_of_exp
        for j=1:number_of_exp
            overlap = min(size(experiments{i}.speed,2),size(experiments{j}.speed,2));
            nmse_speed(i,j) = goodnessOfFit(median(experiments{i}.speed(:,1:overlap),1,'omitnan')', median(experiments{j}.speed(:,1:overlap),1,'omitnan')', 'NMSE');
            nmse_omega(i,j) = goodnessOfFit(median(abs(experiments{i}.omega(:,1:overlap)),1,'omitnan')', median(abs(experiments{j}.omega(:,1:overlap)),1,'omitnan')', 'NMSE');
        end
    end
else
    overlap = min(size(experiments{1}.omega,2),size(experiments{2}.omega,2));
    mse_speed = mean((mean(experiments{1}.speed(:,1:overlap),1,'omitnan')-mean(experiments{2}.speed(:,1:overlap),1,'omitnan')).^2);
    mse_omega = mean((mean(abs(experiments{1}.omega(:,1:overlap)),1,'omitnan')-mean(abs(experiments{2}.omega(:,1:overlap)),1,'omitnan')).^2);
    
    nmse_speed = goodnessOfFit(mean(experiments{2}.speed(:,1:overlap),1,'omitnan')', mean(experiments{1}.speed(:,1:overlap),1,'omitnan')', 'NMSE');
    nmse_omega = goodnessOfFit(mean(abs(experiments{2}.omega(:,1:overlap)),1,'omitnan')', mean(abs(experiments{1}.omega(:,1:overlap)),1,'omitnan')', 'NMSE');
    nmse_total = mean([nmse_speed, nmse_omega]);
    
    nmse_speed_med = goodnessOfFit(median(experiments{2}.speed(:,1:overlap),1,'omitnan')', median(experiments{1}.speed(:,1:overlap),1,'omitnan')', 'NMSE');
    nmse_omega_med = goodnessOfFit(median(abs(experiments{2}.omega(:,1:overlap)),1,'omitnan')', median(abs(experiments{1}.omega(:,1:overlap)),1,'omitnan')', 'NMSE');
    nmse_total_med = mean([nmse_speed_med, nmse_omega_med]);
    
    wmape_speed = mape(median(experiments{2}.speed(:,1:overlap),1,'omitnan'), median(experiments{1}.speed(:,1:overlap),1,'omitnan'),'wMAPE');
    wmape_omega = mape(median(abs(experiments{2}.omega(:,1:overlap)),1,'omitnan'), median(abs(experiments{1}.omega(:,1:overlap)),1,'omitnan'),'wMAPE');
    wmape_total = mean([wmape_speed, wmape_omega]);
        
    % metrics_of_interest = {NMSE_speed, NMSE_omega, NMSE_total}; metrics_tags = ["NMSE_v", "NMSE_\omega", "NMSE_{tot}"];
    % metrics_of_interest = {mape_speed, mape_omega, mape_total}; metrics_tags = ["mape_v", "mape_\omega", "mape_{tot}"];
    metrics_of_interest = {wmape_speed, wmape_omega, wmape_total}; metrics_tags = ["wmape_v", "wmape_\omega", "wmape_{tot}"];
    % metrics_of_interest = {NMSE_total, mape_total, wmape_total}; metrics_tags = ["NMSE_{tot}", "mape_{tot}", "wmape_{tot}"];

    for m=1:length(metrics_tags)
        disp(metrics_tags(m)+" = "+num2str(metrics_of_interest{m},'%.2f'))
    end
end

%% POLTS
if make_plots
    % figure % BOX PLOT - SPEED and ANGULAR VELOCITY - MEAN OVER AGENTS
    % subplot(2,1,1)
    % set(gca,'FontSize',12)
    % for i=1:number_of_exp; data_to_plot{i} = mean(experiments{i}.speed,'omitnan'); end
    % myboxplot(data_to_plot,true, 3, {'b','k','r'})
    % ylabel('$v$ [$\mu$m/s]','Interpreter','Latex','FontSize',16)
    % xticklabels({'REAL','SIMULATED'})
    % title('Average over agents')
    % subplot(2,1,2)
    % set(gca,'FontSize',12)
    % for i=1:number_of_exp; data_to_plot{i} = mean(abs(experiments{i}.omega),'omitnan'); end
    % myboxplot(data_to_plot,true, 3, {'b','k','r'})
    % ylabel('$|\omega|$ [rad/s]','Interpreter','Latex','FontSize',16)
    % xticklabels({'REAL','SIMULATED'})
    % set(gcf,'position',[300,300,300,420])
    % if outputDir
    %     saveas(gcf,fullfile(outputDir, 'comparison_boxplot_meanOnAgents'))
    %     saveas(gcf,fullfile(outputDir, 'comparison_boxplot_meanOnAgents'),'png')
    % end
    
    % figure % BOX PLOT - SPEED and ANGULAR VELOCITY - MEAN OVER TIME
    % subplot(2,1,1)
    % set(gca,'FontSize',12)
    % for i=1:number_of_exp; data_to_plot{i} = mean(experiments{i}.speed,2,'omitnan'); end
    % myboxplot(data_to_plot,true, 3, {'b','k','r'})
    % ylabel('$v$ [$\mu$m/s]','Interpreter','Latex','FontSize',16)
    % xticklabels({'REAL','SIMULATED'})
    % title('Average over time')
    % subplot(2,1,2)
    % set(gca,'FontSize',12)
    % for i=1:number_of_exp; data_to_plot{i} = mean(abs(experiments{i}.omega),2,'omitnan'); end
    % myboxplot(data_to_plot,true, 3, {'b','k','r'})
    % ylabel('$|\omega|$ [rad/s]','Interpreter','Latex','FontSize',16)
    % xticklabels({'REAL','SIMULATED'})
    % set(gcf,'position',[300,300,300,420])
    % if outputDir
    %     saveas(gcf,fullfile(outputDir, 'comparison_boxplot_meanOnTime'))
    %     saveas(gcf,fullfile(outputDir, 'comparison_boxplot_meanOnTime'),'png')
    % end
    
    % figure % BOX PLOT - SPEED and ANGULAR VELOCITY - ALL POINTS
    % subplot(2,1,1)
    % set(gca,'FontSize',12)
    % for i=1:number_of_exp; data_to_plot{i} = experiments{i}.speed(:); end
    % myboxplot(data_to_plot, true, 3, {'b','k','r'})
    % ylabel('$v$ [$\mu$m/s]','Interpreter','Latex','FontSize',16)
    % xticklabels({'REAL','SIMULATED'})
    % title('All points')
    % subplot(2,1,2)
    % set(gca,'FontSize',12)
    % for i=1:number_of_exp; data_to_plot{i} = abs(experiments{i}.omega(:)); end
    % myboxplot(data_to_plot, true, 3, {'b','k','r'})%, [0,0.4470,0.7410])
    % ylabel('$|\omega|$ [rad/s]','Interpreter','Latex','FontSize',16)
    % xticklabels({'REAL','SIMULATED'})
    % set(gcf,'position',[300,300,300,420])
    % if outputDir
    %     saveas(gcf,fullfile(outputDir, 'comparison_boxplot_all'))
    %     saveas(gcf,fullfile(outputDir, 'comparison_boxplot_all'),'png')
    % end
    
    % figure % SCATTER PLOT - SPEED and ANGULAR VELOCITY - MEAN OVER TIME
    % x_to_plot=[];
    % y_to_plot=[];
    % grouping=[];
    % for i=1:number_of_exp; x_to_plot = [x_to_plot,  mean(experiments{i}.speed,2,'omitnan')']; end
    % for i=1:number_of_exp; y_to_plot = [y_to_plot,  mean(abs(experiments{i}.omega),2,'omitnan')']; end
    % for i=1:number_of_exp; grouping = [grouping,  i*ones(1,length(mean(experiments{i}.speed,2,'omitnan')))]; end
    % s=scatterhist(x_to_plot, y_to_plot, 'Location','NorthEast','Direction','out','Group', grouping ,'Color','bk','Kernel','on');%,'Marker','.', 'MarkerSize', 15);
    % xlabel(s,'mean $v$ [$\mu$m/s]','Interpreter','Latex','FontSize',16)
    % ylabel(s,'mean $|\omega|$ [rad/s]','Interpreter','Latex','FontSize',16)
    % s(1).YAxisLocation = 'left';
    % s(1).XAxisLocation = 'bottom';
    % %set(get(gca,'children'),'filled',true)
    % s(2).Position = [0.1    0.82   0.7    0.125];
    % s(3).Position = [0.82   0.1    0.125    0.7];
    % s(1).Position(3) = 0.7;
    % s(1).Position(4) = 0.7;
    % ylim([0,3])
    % xlim([0,120])
    % legend({'REAL','SIMULATED'})
    % if outputDir
    %     saveas(gcf,fullfile(outputDir, 'comparison_scatter_meanOnTime'))
    %     saveas(gcf,fullfile(outputDir, 'comparison_scatter_meanOnTime'),'png')
    % end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%% IN VIVO vs IN SILICO v and w %%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure 
    %%%%%%%%%%%%% Speed plot %%%%%%%%%%%%%%%
    subplot(2,1,1)
    xlim([0,max(timeInstants)])
    ylim([0,120])
    %Plot temporal inputs
    if isvarname('u')
        highlightInputs(timeInstants, u, Render.cmap_inputs(end,:), 0.7);
        highlightInputs(timeInstants, 1-u, Render.cmap_inputs(1,:), 0.7);
    end
    %Plot v in vivo and in silico
%     l1=plot(timeInstants, median(abs(experiments{1}.speed),1,'omitnan'),'linewidth',2,'Color', Render.sim_c);
%     l2=plot(timeInstants, median(abs(experiments{2}.speed),1,'omitnan'),'linewidth',2,'Color', Render.exp_c);
    l2=plotWithShade(timeInstants, median(experiments{2}.speed,1,'omitnan'), quantile(experiments{2}.speed, 0.1, 1), quantile(experiments{2}.speed, 0.9, 1), Render.exp_c, 0.4);
    l1=plotWithShade(timeInstants, median(experiments{1}.speed,1,'omitnan'), quantile(experiments{1}.speed, 0.1, 1), quantile(experiments{1}.speed, 0.9, 1), Render.sim_c, 0.4);
    xlabel('$t$ [s]','Interpreter','Latex','FontSize',16)
    ylabel('$v$ [$\mu$m/s]','Interpreter','Latex','FontSize',16)
    legend('REAL','SIMULATED')
    % legend([l1,l2],'REAL','SIMULATED')
    rng=ylim;
    box on
    %%%%%%%%%%%%% Omega plot %%%%%%%%%%%%%%%
    subplot(2,1,2)
    xlim([0,max(timeInstants)])
    ylim([0,2])
    %Plot the inputs
    if isvarname('u')
        highlightInputs(timeInstants, u, Render.cmap_inputs(end,:), 0.7);
        highlightInputs(timeInstants, 1-u, Render.cmap_inputs(1,:), 0.7);
    end
    %Plot w in vivo and in silico
%     l1=plot(timeInstants, median(abs(experiments{1}.omega),1,'omitnan'),'linewidth',2,'Color', Render.sim_c);
%     l2=plot(timeInstants, median(abs(experiments{2}.omega),1,'omitnan'),'linewidth',2,'Color', Render.exp_c);
    l1=plotWithShade(timeInstants, median(abs(experiments{1}.omega),1,'omitnan'), quantile(abs(experiments{1}.omega), 0.1, 1), quantile(abs(experiments{1}.omega), 0.9, 1), Render.sim_c, 0.3);
    l2=plotWithShade(timeInstants, median(abs(experiments{2}.omega),1,'omitnan'), quantile(abs(experiments{2}.omega), 0.1, 1), quantile(abs(experiments{2}.omega), 0.9, 1), Render.exp_c, 0.3);
    xlabel('$t$ [s]','Interpreter','Latex','FontSize',16)
    ylabel('$|\omega|$ [rad/s]','Interpreter','Latex','FontSize',16)
    % legend([l1,l2],'REAL','SIMULATED')
    rng=ylim;
    box on
    if outputDir
        saveas(gcf,fullfile(outputDir, 'comparison_time_plot'))
        saveas(gcf,fullfile(outputDir, 'comparison_time_plot'),'png')
    end
    

end

end


