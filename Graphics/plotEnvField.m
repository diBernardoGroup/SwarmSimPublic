function [] = plotEnvField(points, values, window, cmap)
    arguments
        points  (1,2)   cell
        values          double
        window  (1,2)   double
        cmap    (:,3)   double = linspace2([1,1,1], [1,0.5,0.5], 100)'
    end
    x_vec = linspace(-window(1)/2,window(1)/2,500);
    y_vec = linspace(-window(2)/2,window(2)/2,500);
    [x_mesh, y_mesh] = meshgrid(x_vec, y_vec);

    F = griddedInterpolant(points,values, 'linear', 'nearest');
    
    imagesc(x_vec,y_vec,F(x_mesh',y_mesh')')
    colormap(cmap)
end

