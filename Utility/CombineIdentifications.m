clear
close all

id_files = ["/Volumes/DOMEPEN/Experiments/comparisons/Euglena_switch_10/combo/identification_OLS_ds3_sign_grad.txt",
            "/Volumes/DOMEPEN/Experiments/comparisons/Euglena_switch_10/combo/identification_OLS_ds1_sign_grad.txt"]; % identification data

output_file = "/Volumes/DOMEPEN/Experiments/comparisons/Euglena_switch_10/combo/identification_OLS_dscombo.txt";

%% LOAD DATA
identifications={};
for i=1:length(id_files)
    identifications{i}=readtable(id_files(i));
end

%% AGGREGATE IDENTIFICATIONS
% combo_id = vertcat(identifications{:});
% writetable(combo_id, output_file ,'Delimiter',' ')

%% COMBINE SPEED AND ANG VEL IDENTIFICATIONS
% use first identification for speed parameters and the second for ang vel
[common_ids, posA, posB] = intersect(identifications{1}.agents,identifications{2}.agents);

combo_id = horzcat(identifications{1}(posA,1:6), identifications{2}(posB,7:end))
writetable(combo_id, output_file ,'Delimiter',' ')

disp(['resulting ',num2str(size(combo_id,1)), ' combo agents from ' ,num2str(size(identifications{1},1)), ' and ', num2str(size(identifications{2},1)), ])

