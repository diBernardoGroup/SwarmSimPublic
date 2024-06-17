function [] = myboxplot(data, significance, whisker, colors)

    arguments
        data
        significance = false
        whisker = 1.5
        colors = get(gca,'ColorOrder')
    end
    
    hold on
    
    %boxplot
    if istable(data)
        boxplot(data.Variables, labels=data.Properties.VariableNames,whisker=whisker,colors='k');
        xlim([0.5 length(data.Properties.VariableNames)+0.5])
    else
        for i=1:length(data)
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
    
    %significance bar
    if significance
        all_data = [];
        for d=1:length(data)
            all_data = [all_data, data{d}'];
        end
        y_max=max(all_data,[],'all');
        y_min=min(all_data,[],'all');
        if length(data)==1
            p = signtest(data{1});
            %[h,p] = ttest(data{1});
        else
            p = ranksum(data{1},data{2});
%             y_max=max(max(data{1}),max(data{2}));
%             y_min=min(min(data{1}),min(data{2}));
            plot([1 1 2 2], [1 1.05 1.05 1]*y_max*1.1, '-k', LineWidth=1)
        end
        
        t=['p=',num2str(p,'%.3f')];
        if p<0.001; t='p<0.001';
        elseif p>0.05; t='ns'; end
        text(mean([1 length(data)]), y_max*1.25, t, HorizontalAlignment='center', FontSize=12)
        ylim([y_min-(y_max-y_min)*0.1, y_max*1.4])
    end
    
end

