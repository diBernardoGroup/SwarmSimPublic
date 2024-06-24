function u = loadInputPattern(data_folder, blur)
    % load input pattern
    if exist(fullfile(data_folder,'patterns_cam/pattern_10.0.jpeg') , 'file')
        inputs=imread(fullfile(data_folder,'patterns_cam/pattern_10.0.jpeg'));
    elseif exist(fullfile(data_folder,'patterns_cam/pattern_30.0.jpeg') , 'file')
        inputs=imread(fullfile(data_folder,'patterns_cam/pattern_30.0.jpeg'));
    else
        inputs=imread(fullfile(data_folder,'patterns_cam/pattern_60.0.jpeg'));
    end
    
    % select blue channel and scale in [0,1]
    %u = flip(double(inputs(:,:,3))',2)/255;    
    u = double(inputs(:,:,3))'/255;    
    
    % blur input pattern
    if blur
        u = movmean2(u,blur);
    end
    
end

