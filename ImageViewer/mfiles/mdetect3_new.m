function fin_ROI=mdetect3_new(man_ROI,inner,outer,image)

man_ROI_margin = bwmorph(man_ROI,'remove'); % for the display of the manually selected area
save_image = image.*(1-man_ROI_margin);
save_image = save_image/max(save_image(:)); % normalization

if ~exist('temporary_saved_variables','dir')
    mkdir('temporary_saved_variables');
end
temp = 'temporary_saved_variables/man_sel_ROI.tiff';
imwrite(save_image,temp,'tif');

outer2 = bwmorph(outer, 'remove');
outer2 = 1 - (outer - outer2);

lower_limit_value = 1.75;
%fin_ROI = man_ROI;
fin_ROI =inner;
fin_area = -1; 
conv_flag = 0;

%figure;

i = 0;

while conv_flag ~= 1
    
    ROI_image = image.*fin_ROI;
    ROI_image_vec = ROI_image(:);
    mean_ROI = mean(ROI_image_vec(ROI_image_vec~=0));
    std_ROI = std(ROI_image_vec(ROI_image_vec~=0));
    lower_limit_ROI = lower_limit_value*std_ROI;
    if exist('thr_old','var')
    %thr_new = (thr_old+mean_ROI - lower_limit_ROI)/2;
       thr_new = mean_ROI - lower_limit_ROI;
    else
        thr_new = mean_ROI - lower_limit_ROI;
    end
    BW2 = roicolor(image,thr_new, max(max(image))); % select regions of interest based on lower cut-off criterion
    
    thr_old =thr_new;
    
    labeled_mask = bwlabel(BW2); % labeling of different created masks
    ROI_label = labeled_mask.*fin_ROI; % only the labeled areas in the inner area

    Y = sort(ROI_label(:)); % sorted ROI_label
    N = histc(ROI_label(:),Y); % contains the "areas" of the labeled areas
    Y(N==0) = [];
    N(N==0) = [];
    N(Y==0)=0;
    [largest_area,index] = max(N); % search for the largest area
    label = Y(index);
    ROI_new = zeros(size(image,1), size(image,2));
    ROI_new(labeled_mask == label) = 1;
    %ROI_new = imfill(ROI_new);
    %fin_ROI = (ROI_new | inner) & outer; 
    fin_ROI =ROI_new;
    area_diff = fin_area-sum(fin_ROI(:));
    fin_area = sum(fin_ROI(:));

    if area_diff == 0
         conv_flag = 1;
    end
     i=i+1;
     
     
     
    if i<25
        if conv_flag == 1
          %  subplot(5,5,i); imshow(fin_ROI + outer2-bwmorph(inner, 'remove'),[]); 
           % title('Convergence')
           % subplot(5,5,i+1); imshow(fin_ROI.*image ,[]); title('The final ROI')
           % savename = ['temporary_saved_variables/' 'ROI_algorithm_' num2str(100*lower_limit_value) '_std.png']; % save as fig...
           % saveas(gcf,savename)
        else
           % subplot(5,5,i); imshow(fin_ROI + outer2-bwmorph(inner, 'remove'),[]); 
        end
    elseif i<50
        if ~exist('f_handle2','var')
        % f_handle2 = figure;
        end
        
        if conv_flag == 1
           %figure;
           %subplot(5,5,i-24); imshow(fin_ROI + outer2-bwmorph(inner, 'remove'),[]); %title([num2str(mean_ROI) '-' num2str(std_ROI)]);
           %title('Convergence')
           % subplot(5,5,i-24+1); imshow(fin_ROI.*image ,[]); title('The final ROI')
           % savename = ['temporary_saved_variables/' 'ROI_algorithm_' num2str(100*lower_limit_value) '_std2.png']; % save as fig...
           % saveas(gcf,savename)
        else
           % subplot(5,5,i-24); imshow(fin_ROI + outer2-bwmorph(inner, 'remove'),[]); %title([num2str(mean_ROI) '-' num2str(std_ROI)]);
        end
        
    else
        disp([mfilename,': ROI detection failed to converge within 50 iterations']);
        return;
    end
   
end

if sum(fin_ROI(:)) - sum(outer(:)) == 0
         warndlg('Dectected roi is the same as the outer boundary'); 
end
