function [] = myboxplot(data)
    hold on
    
    %boxplot
    for i=1:length(data)
        boxplot(data{i},positions=i,labels=num2str(i))
    end
    
    %significance bar
    p = ranksum(data{1},data{2});
    y_max=max(max(data{1}),max(data{2}));
    y_min=min(min(data{1}),min(data{2}));
    plot([1 1 2 2], [1 1.05 1.05 1]*y_max*1.1, '-k')
    t=['p=',num2str(p,'%.3f')];
    if p<0.001; t='p<0.001'; 
    elseif p>0.05; t='ns'; end
    text(mean([1 2]), y_max*1.25, t, HorizontalAlignment='center')
    
    xticks(1:length(data))
    xlim([0.5 length(data)+0.5])
    ylim([y_min*0.8, y_max*1.4])
end

