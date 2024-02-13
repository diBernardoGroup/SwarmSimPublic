



%data_folder = '/Volumes/DOMEPEN/Experiments/2023_06_15_Euglena_1/tracking_2023_10_12';  % off
%data_folder = '/Volumes/DOMEPEN/Experiments/2023_06_15_Euglena_7/tracking_2023_10_16';  % switch10s
%data_folder = '/Volumes/DOMEPEN/Experiments/2023_06_26_Euglena_37/tracking_2023_10_12'; % circle light
%data_folder = '/Volumes/DOMEPEN/Experiments/2023_07_10_Euglena_26/tracking_2024_01_30'; % circle light high denisty
%data_folder = '/Volumes/DOMEPEN/Experiments/2023_07_10_Euglena_21/tracking_2024_01_30'; % circle dark
%data_folder = '/Volumes/DOMEPEN/Experiments/2023_06_26_Euglena_33/tracking_2023_10_12'; % gradient central light
data_folder = '/Volumes/DOMEPEN/Experiments/2023_06_26_Euglena_34/tracking_2023_10_12'; % gradient central dark

id_folder = '/Volumes/DOMEPEN/Experiments/2023_06_15_Euglena_7/tracking_2023_10_16'; % folder with identification data


identification_file_names = ["identification_OLS_ds1_sign.txt","identification_OLS_ds2_sign.txt","identification_OLS_ds3_sign.txt";
                             %"identification_OLS_ds1_abs.txt","identification_OLS_ds2_abs.txt","identification_OLS_ds3_abs.txt";
                             "identification_GB_ds1_sign.txt","identification_GB_ds2_sign.txt","identification_GB_ds3_sign.txt";
                             "identification_GBCT_ds1_sign_grad.txt","identification_GBCT_ds2_sign_grad.txt","identification_GBCT_ds3_sign_grad.txt";
                             "identification_GBCT_ds1_sign_diff.txt","identification_GBCT_ds2_sign_diff.txt","identification_GBCT_ds3_sign_diff.txt"];

tags = ["OLS","GBDT","GBCT grad","GBCT diff"];
                         
identifications={};
for i=1:size(identification_file_names,2)   % for each down sampling value
for j=1:size(identification_file_names,1)   % for each technique
    identifications{j,i}=readtable(fullfile(id_folder,identification_file_names(j,i)));
end
end


%% PRINT RESULTS

fprintf(join(['Tech','DS',repmat("%s",1,length(identifications{1,1}.Properties.VariableNames)-5),'\n'],'\t'),identifications{j,1}.Properties.VariableNames{2:end-4})
for j=1:size(identification_file_names,1) % for each technique
    disp(tags(j))
    for i=1:size(identification_file_names,2) % for each down sampling value
        fprintf(join(['',num2str(i),repmat("%.2f",1,length(identifications{j,i}.Properties.VariableNames)-5),'\n'],'\t'),mean(identifications{j,i}{:,2:end-4}))
    end
end

figure
colors = get(gca,'ColorOrder');
for k=1:10 % for each parameter
    ax=subplot(2,5,k);
    
    for j=1:size(identification_file_names,1) % for each technique
        medians = [median(identifications{j,1}{:,k+1}), median(identifications{j,2}{:,k+1}), median(identifications{j,3}{:,k+1})];
        quartiles1 = [quantile(identifications{j,1}{:,k+1},0.25), quantile(identifications{j,2}{:,k+1},0.25), quantile(identifications{j,3}{:,k+1},0.25)];
        quartiles3 = [quantile(identifications{j,1}{:,k+1},0.75), quantile(identifications{j,2}{:,k+1},0.75), quantile(identifications{j,3}{:,k+1},0.75)];
        line=plotWithShade([1,2,3],medians,quartiles1,quartiles3, colors(j,:), 0.3);
    end
    legend({'',tags(1),'',tags(2),'',tags(3),'',tags(4)})
    xticks([1,2,3])
    title(identifications{1}.Properties.VariableNames(k+1))
end

