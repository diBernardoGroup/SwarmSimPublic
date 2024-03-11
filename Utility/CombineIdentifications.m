clear
close all

id_files = ["/Volumes/DOMEPEN/Experiments/2023_06_15_Euglena_7/tracking_2023_10_16/identification_OLS_ds1_sign.txt",
            "/Volumes/DOMEPEN/Experiments/2023_06_15_Euglena_7/tracking_2023_10_16/identification_OLS_ds2_sign.txt"]; % identification data

output_file = "/Volumes/DOMEPEN/Experiments/2023_06_15_Euglena_7/tracking_2023_10_16/identification_combo.txt";

%% LOAD DATA
identifications={};
for i=1:length(id_files)
    identifications{i}=readtable(id_files(i));
end

combo_id = vertcat(identifications{:})

writetable(combo_id, output_file ,'Delimiter',' ')


