function u = loadInputPattern(data_folder, blur, timeInstant)

arguments
    data_folder     char
    blur            double 
    timeInstant     double = 180
end

    % get patterns in the input folder
    pattern_names = strsplit(ls(fullfile(data_folder,'patterns_cam')),{'\t','\n'});
    pattern_names = pattern_names(1:end-1);
    
    for i=1:length(pattern_names)
        pattern_instants(i) = sscanf(pattern_names{i},'pattern_%f.jpeg');
    end
    
    % last pattern before timeInstant
    my_instant = max(pattern_instants(pattern_instants<=timeInstant));
    my_pattern_name = sprintf('pattern_%04.1f.jpeg',my_instant);
    
    % load pattern
    inputs=imread(fullfile(data_folder,'patterns_cam',my_pattern_name));
    
    % select blue channel and scale in [0,1]
    %u = flip(double(inputs(:,:,3))',2)/255;    
    u = double(inputs(:,:,3))'/255;    
    
    % blur input pattern
    if blur
        u = movmean2(u,blur);
    end
    
end

