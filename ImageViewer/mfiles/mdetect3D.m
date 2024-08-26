function roi = mdetect3D(roi,inner,outer,image)
% Now this function can deal with three dimensional data as well.
% man_ROI,inner,outer,image are matrices of the same size, containing data
% for roi, inner bound, outer bound and image, repectively.

lower_limit_value = 1.75;

%% Finding ROI: Take the largest labeled area inside the manually selected area

fin_area = -1; 
conv_flag = 0;
i = 0;
h = waitbar(0,'3d ROI detection');
while conv_flag ~= 1
    
    waitbar(i/50,h);
    ROI_image = image.*roi;
    ROI_image_vec = ROI_image(:);
    mean_roi = mean(ROI_image_vec(ROI_image_vec~=0));
    std_roi = std(ROI_image_vec(ROI_image_vec~=0));
    lower_limit_roi = lower_limit_value*std_roi;

    BW2 = (image>(mean_roi - lower_limit_roi)); % select regions of interest based on lower cut-off criterion
    labeled_mask = bwlabeln(BW2); % labeling of different created masks
    ROI_label = labeled_mask.*roi;

    
    Y = sort(ROI_label(:)); % sorted ROI_label
    N = histc(ROI_label(:),Y); % contains the "areas" of the labeled areas
    Y(N==0) = [];
    N(N==0) = [];
    N(Y==0)=0;
    [largest_area,index] = max(N); % search for the largest area
    label = Y(index);
    ROI_new = zeros(size(image));
    ROI_new(labeled_mask == label) = 1;
    for ifl = 1:size(roi,3)
     ROI_new(:,:,ifl) = imfill(ROI_new(:,:,ifl));
    end
    roi = (ROI_new | inner) & outer; 
    
    area_diff = fin_area-sum(roi(:));
    fin_area = sum(roi(:));

    if area_diff == 0
        conv_flag = 1;
    end
    
    i=i+1;
    
    if i>=50 && conv_flag == 0
        warndlg([mfilename,': ROI detection failed to converge within 50 iterations']);
        delete(h);
        return;
    end
    
end

delete(h);
