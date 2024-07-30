function [] = myboxplot(data, significance, whisker, colors, out_threshold)

arguments
    data
    significance = false
    whisker = 1.5
    colors = get(gca,'ColorOrder')
    out_threshold = 2
end

hold on

%boxplot
if istable(data)
    boxplot(data.Variables, labels=data.Properties.VariableNames,whisker=whisker,colors='k');
    xlim([0.5 length(data.Properties.VariableNames)+0.5])
else
    for i=1:length(data)
        data{i}=data{i}(~isnan(data{i}));
        out_indx = isoutlier(data{i},'quartiles',thresholdfactor=out_threshold);
        data{i}=data{i}(~out_indx);
        boxplot(data{i},positions=i,labels=num2str(i),whisker=whisker,colors='k');
    end
    xticks(1:length(data))
    xlim([0.5 length(data)+0.5])
end

h = findobj(gcf,'tag','Outliers');
for j=1:length(h)
    h(j).MarkerEdgeColor = validatecolor('w')*0.5;
end

h = findobj('LineStyle','--'); set(h, 'LineStyle','-');
h = findobj(gca,'Type','line'); set(h, 'LineWidth',1);
h = flip(findobj(gca,'Tag','Box'));
if iscell(colors) && ischar(colors{1})
    colors=validatecolor(colors,'multiple');
end
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(mod(j-1,7)+1,:),'FaceAlpha',.5);
end

% stat test and significance bar
if significance
    all_data = [];
    for d=1:length(data)
        if size(data{d},1)>size(data{d},2); data{d}=data{d}'; end
        all_data = [all_data, data{d}];
    end
    y_max=max(all_data,[],'all');
    y_min=min(all_data,[],'all');
    if y_min==0 && y_max==0
        y_max=1;
        y_min=-1;
    end
    
    
    if length(data)==1      % single distribution test
        p = signtest(data{1});
        %[h,p] = ttest(data{1});
        t=['p=',num2str(p,'%.3f')];
        if p<0.001; t='p<0.001';
        elseif p>0.05; t='ns'; end
        text(d+0.5, y_max*1.25, t, HorizontalAlignment='center', FontSize=12)
        
    else                    % two distributions test
        for d=1:length(data)-1
            % p = ranksum(data{d},data{d+1});
            [~,p]=kstest2(data{d},data{d+1});
            plot([0 0 1 1]+d, [1 1.04 1.04 1]*y_max+(y_max-y_min)*0.025, '-k', LineWidth=1)
            
            t=['p=',num2str(p,'%.3f')];
            if p<0.001; t='p<0.001';
            elseif p>0.05; t='ns'; end
            text(d+0.5, y_max+(y_max-y_min)*0.11, t, HorizontalAlignment='center', FontSize=12)
        end
    end
    
    % automatic yaxis scaling
    ylim([y_min-(y_max-y_min)*0.1, y_max+(y_max-y_min)*0.15])
end

end

