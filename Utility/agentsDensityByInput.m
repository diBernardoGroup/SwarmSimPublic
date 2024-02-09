function [density_by_input, bins, norm_slope, c_coeff, coefficents] = agentsDensityByInput(points, values, x, window)
    F = griddedInterpolant(points,values, 'linear', 'nearest');
    
    if size(x,2)>2
        mask=x;
        x_vec = linspace(window(1),window(2),size(mask,2));
        y_vec = linspace(window(3),window(4),size(mask,1));
        [x_mesh, y_mesh] = meshgrid(x_vec, y_vec);
        envInput = F(x_mesh(mask),y_mesh(mask));
    else
        x_vec = linspace(window(1),window(2),100);
        y_vec = linspace(window(3),window(4),100);
        [x_mesh, y_mesh] = meshgrid(x_vec, y_vec);
        envInput = F(x(:,1),x(:,2));    % input intensity measured by the agents
    end
    
    [pixels_by_input,bins] = histcounts(F(x_mesh',y_mesh'), 3);
    [agents_by_input,bins] = histcounts(envInput, bins);
    density_by_input = agents_by_input./pixels_by_input;
    c_matrix = corrcoef(density_by_input, bins(1:end-1));
    c_coeff = c_matrix(1,2);
    
    coefficents = [ones(length(bins(1:end-1)),1),(bins(1:end-1)+bins(2:end))'/2]\density_by_input';
    norm_slope = coefficents(2)/abs(mean(density_by_input));
end

