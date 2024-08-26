function fin_ROI = adjustROI(fin_ROI,tissue_BW,c_image)

    temp = c_image.*tissue_BW;
    mean_tissue_uptake = mean(temp(:));
    std_tissue_uptake = std(temp(:));
    uptake_threshold = mean_tissue_uptake + 2*std_tissue_uptake;
    
    uptake_mask=(c_image>=uptake_threshold);
    fin_ROI = fin_ROI.*uptake_mask;