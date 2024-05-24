clear
close all

% id_folder = '/Volumes/DOMEPEN/Experiments/2023_06_15_Euglena_7/tracking_2023_10_16'; % folder with identification data
id_folder = '/Volumes/DOMEPEN/Experiments/comparisons/Euglena_switch_10/combo5';

% identification_file_names = ["identification_OLS_ds1_sign.txt","identification_OLS_ds2_sign.txt","identification_OLS_ds3_sign.txt";
%                              %"identification_OLS_ds1_abs.txt","identification_OLS_ds2_abs.txt","identification_OLS_ds3_abs.txt";
%                              "identification_GB_ds1_sign.txt","identification_GB_ds2_sign.txt","identification_GB_ds3_sign.txt";
%                              "identification_GBCT_ds1_sign_grad.txt","identification_GBCT_ds2_sign_grad.txt","identification_GBCT_ds3_sign_grad.txt";
%                              "identification_GBCT_ds1_sign_diff.txt","identification_GBCT_ds2_sign_diff.txt","identification_GBCT_ds3_sign_diff.txt"];
                         
identification_file_names = ["identification_OLS_ds1_sign_grad.txt","identification_OLS_ds2_sign_grad.txt","identification_OLS_ds3_sign_grad.txt";
                             "identification_GBCT_ds1_sign_grad.txt","identification_GBCT_ds2_sign_grad.txt","identification_GBCT_ds3_sign_grad.txt"];
tags = ["OLS","GBDT","GBCT grad","GBCT diff"];

identification_file_names = ["identification_OLS_ds1.txt","identification_OLS_dscombo.txt","identification_OLS_ds3.txt"];

identification_file_names = ["identification_OLS+GB_ds1_old.txt"; 
                                "identification_OLS+GB_ds1_diff.txt"; 
                                "identification_OLS+GB_ds1_diff_min20.txt";
                                "identification_OLS+GB_ds1_diff_noalpha.txt";
                                "identification_OLS+GB_ds1_diff_nosign.txt";
                                "identification_OLS+GB_ds1_diff_median";
                                "identification_OLS+GB_ds1_diff_smooth_median";
                                "identification_OLS+GB_ds1_diff_median_smooth";
                                "identification_OLS+GB_ds1_diff_smooth_median_smooth"];
tags = ["old","min=5s","min=20s","no\_alpha","nosign","med","smooth_med","med_smooth","s_med_s"];

identification_file_names =    ["identification_OLS+GB_ds1_diff_sign.txt"; 
                                "identification_OLS+GB_ds1_diff_nosign.txt";
                                "identification_OLS+GB_ds1_diff_sign_nomu.txt";
                                "identification_OLS+GB_ds1_diff_nosign_nomu.txt"];
tags = ["signed","nosign","signed-nomu","nosign-nomu"];
ref_sim = [];

dT = 0.01;
deltaT = 0.5;

%% LOAD DATA
identifications={};
for i=1:size(identification_file_names,2)   % for each down sampling value
for j=1:size(identification_file_names,1)   % for each technique
    identifications{j,i}=readtable(fullfile(id_folder,identification_file_names(j,i)));
end
end

speed=load(fullfile(id_folder,'speeds_smooth.txt'));
omega=load(fullfile(id_folder,'ang_vel_smooth.txt'));
inputs=load(fullfile(id_folder,'inputs.txt'));
timeInstants = [0:size(speed,1)-1] * deltaT;
u=inputs(:,1)/255;
u_dot = [0;diff(u)]/deltaT;
%u_dot = gradient(u)/deltaT;
u_dot = max(u_dot,0);

%% simulate average behaviour
t_sim=0:dT:max(timeInstants);
s_sim=nan(length(t_sim),1);
w_sim=nan(length(t_sim),1);
u_sim=nan(length(t_sim),2);
for j=1:size(identification_file_names,1)       % for each technique
    for i=1:size(identification_file_names,2)   % for each down sampling value
        theta_s_mean=median(identifications{j,i}.theta_s);
        mu_s_mean   =median(identifications{j,i}.mu_s);
        alpha_s_mean=median(identifications{j,i}.alpha_s);
        beta_s_mean =median(identifications{j,i}.beta_s);
        theta_w_mean=median(identifications{j,i}.theta_w);
        mu_w_mean   =median(abs(omega),'all','omitnan');
        alpha_w_mean=median(identifications{j,i}.alpha_w);
        beta_w_mean =median(identifications{j,i}.beta_w);
        s_sim(1)=mu_s_mean;
        w_sim(1)=mu_w_mean;
        for t=1:length(t_sim)-1                 % for each time instant
            u_sim(t,:)= [u(ceil(t*dT/deltaT)),u_dot(ceil(t*dT/deltaT))];
            s_sim(t+1)= s_sim(t) + (theta_s_mean * (mu_s_mean-s_sim(t)) + alpha_s_mean * u_sim(t,1) + beta_s_mean * u_sim(t,2) ) *dT;
            w_sim(t+1)= w_sim(t) + (theta_w_mean * (mu_w_mean-w_sim(t)) + alpha_w_mean * u_sim(t,1) + beta_w_mean * u_sim(t,2) ) *dT;
        end
        
        nmse_med_speed(j,i) = goodnessOfFit(interp1(t_sim,s_sim,timeInstants)', median(speed,2,'omitnan'), 'NMSE');
        nmse_med_omega(j,i) = goodnessOfFit(interp1(t_sim,w_sim,timeInstants)', median(abs(omega),2,'omitnan'), 'NMSE');
        nmse_med_total(j,i) = mean([nmse_med_speed(j,i), nmse_med_omega(j,i)]);
        
        wmape_speed(j,i) = mape(interp1(t_sim,s_sim,timeInstants)', median(speed,2,'omitnan'),'wMAPE');
        wmape_omega(j,i) = mape(interp1(t_sim,w_sim,timeInstants)', median(abs(omega),2,'omitnan'),'wMAPE');
        wmape_total(j,i) = mean([wmape_speed(j,i), wmape_omega(j,i)]);
    end
end



%% PRINT RESULTS

% metrics_of_interest = {nmse_speed, nmse_omega, nmse_total}; metrics_tags = ["nmse_v", "nmse_\omega", "nmse_{tot}"];
% metrics_of_interest = {mape_speed, mape_omega, mape_total}; metrics_tags = ["mape_v", "mape_\omega", "mape_{tot}"];
metrics_of_interest = {wmape_speed, wmape_omega, wmape_total}; metrics_tags = ["wmape_v", "wmape_\omega", "wmape_{tot}"];
% metrics_of_interest = {nmse_total, mape_total, wmape_total}; metrics_tags = ["nmse_{tot}", "mape_{tot}", "wmape_{tot}"];
metrics_color = ['b','r','k'];

fprintf(join(['Tech','DS',repmat("%s",1,length(identifications{1,1}.Properties.VariableNames)-5)],'\t'),identifications{j,1}.Properties.VariableNames{2:end-4})
fprintf(join(['\t',repmat("%s",1,length(metrics_tags)),'\n'],'\t'),metrics_tags)
for j=1:size(identification_file_names,1) % for each technique
    disp(tags(j))
    for i=1:size(identification_file_names,2) % for each down sampling value
        % print average parameters value
        fprintf(join(['',num2str(i),repmat("%.2f",1,length(identifications{j,i}.Properties.VariableNames)-5)],'\t'),mean(identifications{j,i}{:,2:end-4}))
        % print metrics
        for k=1:length(metrics_of_interest); data_to_print(k) = metrics_of_interest{k}(j,i); end
        fprintf(join(['',repmat("%.2f",1,length(metrics_of_interest))],'\t\t'),data_to_print)
        fprintf('\n')
    end
end

if size(identification_file_names,2)>1 % if considering multiple down sampling values
figure % Parameters
colors = get(gca,'ColorOrder');
for k=1:10 % for each parameter
    ax=subplot(2,5,k);
    
    for j=1:size(identification_file_names,1) % for each technique
        medians = [median(identifications{j,1}{:,k+1}), median(identifications{j,2}{:,k+1}), median(identifications{j,3}{:,k+1})];
        quartiles1 = [quantile(identifications{j,1}{:,k+1},0.25), quantile(identifications{j,2}{:,k+1},0.25), quantile(identifications{j,3}{:,k+1},0.25)];
        quartiles3 = [quantile(identifications{j,1}{:,k+1},0.75), quantile(identifications{j,2}{:,k+1},0.75), quantile(identifications{j,3}{:,k+1},0.75)];
        line=plotWithShade([1,2,3],medians,quartiles1,quartiles3, colors(j,:), 0.3);
    end
    legend({'',tags(1),'',tags(2)})
    xticks([1,2,3])
    title(identifications{1}.Properties.VariableNames(k+1))
end

figure % NMSE
rng = max([nmse_med_speed,nmse_mean_speed,nmse_med_omega,nmse_mean_omega],[],'all');
for j=1:size(identification_file_names,1) % for each technique
    subplot(3,2,1); hold on
    plot([1,2,3], nmse_mean_speed(j,:), 'color', colors(j,:), 'marker', 'o')
    title('NMSE mean s')
    legend(tags(1:j))
    xticks([1,2,3])
    ylim([0 rng])
    subplot(3,2,2); hold on
    plot([1,2,3], nmse_med_speed(j,:), 'color', colors(j,:), 'marker', 'o')
    title('NMSE med s')
    legend(tags(1:j))
    xticks([1,2,3])
    ylim([0 rng])
   
    subplot(3,2,3); hold on
    plot([1,2,3], nmse_mean_omega(j,:), 'color', colors(j,:), 'marker', 'o')
    title('NMSE mean w')
    legend(tags(1:j))
    xticks([1,2,3])
    ylim([0 rng])
    subplot(3,2,4); hold on
    plot([1,2,3], nmse_med_omega(j,:), 'color', colors(j,:), 'marker', 'o')
    title('NMSE med w')
    legend(tags(1:j))
    xticks([1,2,3])
    ylim([0 rng])
    
    subplot(3,2,5); hold on
    plot([1,2,3], nmse_mean_total(j,:), 'color', colors(j,:), 'marker', 'o')
    title('NMSE mean tot')
    legend(tags(1:j))
    xticks([1,2,3])
    ylim([0 rng])
    subplot(3,2,6); hold on
    plot([1,2,3], nmse_med_total(j,:), 'color', colors(j,:), 'marker', 'o')
    title('NMSE med tot')
    legend(tags(1:j))
    xticks([1,2,3])
    ylim([0 rng])
end

else % if considering a single down sampling value
    
        figure('Position',[100 100 1900 600]); % PARAMETERS BOXPLOTS
    for k=1:10 % for each parameter
        ax=subplot(2,5,k);
        
        for j=1:size(identification_file_names,1); data_to_plot{j} = identifications{j,1}{:,k+1}; end
        myboxplot(data_to_plot, false, 3);
        xticks([])
        set(ax,'PositionConstraint','innerposition')
        yline(0,'Color',[0.5,0.5,0.5])
        %l=max(identification.(i+1))*1.1;ylim([-l/15,l]);yticks([0:l/3:l])
        xlabel(identifications{1}.Properties.VariableNames(k+1))
        set(gca,'FontSize',16)
    end
    legend(tags)
    
    figure('Position',[100 100 1900 600]); % metrics scatter comparison
    hold on
    for k=1:length(metrics_of_interest)
        plots(k,:)=scatter([1:length(tags)]-(length(metrics_of_interest)-1)*0.1+(k-1)*0.2,metrics_of_interest{k},50,metrics_color(k));
        scatter([1:length(tags)]-(length(metrics_of_interest)-1)*0.1+(k-1)*0.2,metrics_of_interest{k}(:,1),50,metrics_color(k),"filled");
    end
    xticks([1:length(tags)])
    xticklabels(tags)
    set(gca, 'TickLabelInterpreter', 'none');
    xlim([0,length(tags)+1])
    ylim([0, max([metrics_of_interest{:}],[],'all')*1.1])
    legend(plots(:,1),metrics_tags)
    set(gca,'FontSize',14)
    box on
    set(gca,'XGrid','off','YGrid','on')
    
end

