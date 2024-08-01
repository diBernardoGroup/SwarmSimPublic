function mask = detectObjects(data_folder, background_sub, brightness_thresh,timeInstant, brightness_scaling)


arguments
    data_folder           char
    background_sub        logical
    brightness_thresh     double 
    timeInstant           double        = 180
    brightness_scaling    logical       = true
end



    % load last image
    if timeInstant<10
        img=imread(fullfile(data_folder,sprintf('images/fig_0%.1f.jpeg',timeInstant)));
    else
        img=imread(fullfile(data_folder,sprintf('images/fig_%.1f.jpeg',timeInstant)));
    end
    
    img_grey = double(img(:,:,1))/255;
    
    % adjust brightness to use the full dynamic range
    if brightness_scaling
    min_val = min(img_grey(360:720,640:1280),[],'all');
    img_grey = max(min_val, img_grey);
    img_grey = img_grey - min_val;
    max_val = max(img_grey(360:720,640:1280),[],'all');
    img_grey = img_grey/max_val;
    end
    
    % background subtraction
    if background_sub
        background = double(imread(fullfile(data_folder,'background.jpeg')))/255;
        foreground = max(img_grey-background,0);
    else
        foreground = img_grey;
    end

    % theresholding
    mask=flip(foreground>brightness_thresh);

end

