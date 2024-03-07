function [c_coeff,norm_slope,coefficents] = linearDependence(x,y)

    nanIdx = isnan(x) | isnan(y);
    x=x(~nanIdx);
    y=y(~nanIdx);
    
    c_matrix = corrcoef(x, y, 'rows', 'complete');
    c_coeff = c_matrix(1,2);
    
    coefficents = [ones(length(x),1), x]\y;
    norm_slope = coefficents(2)/abs(mean(y))*abs(mean(x));
    
%     figure;scatter(x,y);hold on;plot([0,max(x)], coefficents(1)+coefficents(2)*[0,max(x)]); yline(0,color=[1,1,1]*0.75)
%     text(max(x),max(y)*1,['\rho=',num2str(c_coeff,3)],'HorizontalAlignment','right','FontSize',14)
%     text(max(x),max(y)*0.9,['norm slope=',num2str(norm_slope,3)],'HorizontalAlignment','right','FontSize',14)
%     box
end

