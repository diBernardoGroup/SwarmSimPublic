function line= plotWithShade(x,yLine,yShadeBottom,yShadeTop, color, alpha)
%
%plotWithShade polts a solid line and a shaded area.
%   Used to plot the average and the min-max range of a value.
%   Point markers are added if the line has no more then 30 points.
%
%   line= plotWithShade(x,yLine,yShadeBottom,yShadeTop, color, alpha)
%
%   Inputs:
%       x are the x values (vector)
%       yLine are the y values of the solid line (vector same length of x)
%       yShadeBottom are the y values of the bottom of the shaded area (vector)
%       yShadeTop are the y values of the top of the shaded area (vector)
%       color is the color of the line and the shaded area (color string)
%       alpha is the opacity of the shaded area (scalar)
%
%   Outputs:
%       line is the plot of the solid line
%
%   Authors:    Andrea Giusti and Gian Carlo Maffettone
%   Date:       2022
%

    yLine=yLine(~isnan(yLine));
    yShadeBottom=yShadeBottom(~isnan(yLine));
    yShadeTop=yShadeTop(~isnan(yLine));
    x=x(~isnan(yLine));
    
    number_of_points=length(x);
    
    if size(yShadeBottom,1) == 1
        yShadeBottom=yShadeBottom';
        yShadeTop=yShadeTop';
    end
    
    assert(length(yLine)==number_of_points, "yLine must have the same length of x")
    assert(length(yShadeBottom)==number_of_points, "yShadeBottom must have the same length of x")
    assert(length(yShadeTop)==number_of_points, "yShadeTop must have the same length of x")
    
    if number_of_points>30
       mark= 'none';
    else
       mark= 'o';
    end
    
    hold on
    patch([x fliplr(x)], [yShadeBottom; flipud(yShadeTop)]', color, 'FaceAlpha',alpha, 'LineStyle','none')
    line=plot(x, yLine, 'Marker',mark,'Color',color,'LineWidth',2,'MarkerSize',6);
end

