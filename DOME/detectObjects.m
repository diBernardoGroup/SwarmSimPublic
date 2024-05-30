function mask = detectObjects(data_folder, background_sub, brightness_thresh)

    % load last image
    img=imread(fullfile(data_folder,'images/fig_180.0.jpeg'));
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

