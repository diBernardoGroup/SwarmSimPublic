function mask = detectObjects(data_folder, background_sub, brightness_thresh,timeInstant)


arguments
    data_folder           char
    background_sub        logical
    brightness_thresh     double 
    timeInstant           double        = 180
end



    % load last image
    if timeInstant<10
        img=imread(fullfile(data_folder,sprintf('images/fig_0%.1f.jpeg',timeInstant)));
    else
        img=imread(fullfile(data_folder,sprintf('images/fig_%.1f.jpeg',timeInstant)));
    end
    
    img_grey = double(img(:,:,1))/255;

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

