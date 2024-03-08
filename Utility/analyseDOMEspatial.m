function [mask, u] = analyseDOMEspatial(data_folder, background_sub, brightness_thresh)

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

    % load input pattern
    if exist(fullfile(data_folder,'patterns_cam/pattern_10.0.jpeg') , 'file')
        inputs=imread(fullfile(data_folder,'patterns_cam/pattern_10.0.jpeg'));
    else
        inputs=imread(fullfile(data_folder,'patterns_cam/pattern_30.0.jpeg'));
    end
    u = flip(double(inputs(:,:,3))',2)/255;    %select blue channel and scale in [0,1]
    u = movmean2(u,51);

end

