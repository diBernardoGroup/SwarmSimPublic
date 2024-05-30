function [] = myHistogram(data,edges,kest)

%Plots an histogram given the edges and gives you the possibility to plot
%an estimated pdf using Kernel Density Estimation

arguments
    data
    edges=[];
    kest=false;
end

if kest == false

    if isempty(edges)
        histogram(data,'Normalization','pdf');
    else
        histogram(data,edges,'Normalization','pdf');
    end

else

    [F, XF] = ksdensity(data);
    plot(XF,F,'LineWidth',2);

end


end