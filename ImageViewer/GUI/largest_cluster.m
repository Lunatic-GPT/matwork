function roi=largest_cluster(mask)   

    labeled_mask = bwlabel(mask>0); % labeling of different created masks
    ROI_label = labeled_mask;

    Y = sort(ROI_label(:)); % sorted ROI_label
    N = histc(ROI_label(:),Y); % contains the "areas" of the labeled areas
    Y(N==0) = [];
    N(N==0) = [];
    N(Y==0)=0;
    [largest_area,index] = max(N); % search for the largest area
    
    label = Y(index);
    roi = zeros(size(mask));
    roi(labeled_mask == label) = 1;
    