function [density_by_input, bins, norm_slope, c_coeff, coefficents, agents_by_input,pixels_by_input, u_values] = agentsDensityByInput(points, values, x, window)
    F = griddedInterpolant(points,values, 'linear', 'nearest');
    
    if size(x,2)>2
        mask=x;
        x_vec = linspace(window(1),window(2),size(mask,2));
        y_vec = linspace(window(3),window(4),size(mask,1));
        [x_mesh, y_mesh] = meshgrid(x_vec, y_vec);
        u_values = F(x_mesh(mask),y_mesh(mask));
    else
        x_vec = linspace(window(1),window(2),1000);
        y_vec = linspace(window(3),window(4),1000);
        [x_mesh, y_mesh] = meshgrid(x_vec, y_vec);
        u_values = F(x(:,1),x(:,2));    % input intensity measured by the agents
    end
    
    [pixels_by_input,bins] = histcounts(F(x_mesh',y_mesh'), [0:1/5:1]);
    [agents_by_input,bins] = histcounts(u_values, bins);
    density_by_input = agents_by_input./pixels_by_input;
    density_by_input = density_by_input/sum(density_by_input);
    [c_coeff,norm_slope,coefficents] = linearDependence((bins(1:end-1)+bins(2:end))'/2, density_by_input');
end

