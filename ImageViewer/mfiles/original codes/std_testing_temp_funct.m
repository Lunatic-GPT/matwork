function std_testing_temp_funct
%% Manually created mask 
load Research_workspace
load 'C:\Documents and Settings\hahntobi\My Documents\MATLAB\temporary saved variables\BW_cut' BW_cut;

image = image_out(2).*(1-BW_cut);

%roi = [162.2 296.5; 150.3 286.9;  156.2 276.2; 167 283.4];
roi = pos;

man_ROI = roipoly(image, roi(:,1)', roi(:,2)'); %create mask
% figure(2);
% imshow(man_ROI);
%temp_im = image_out(1).*uint8(man_ROI); % crop image to mask area
temp_im = image.*man_ROI;





%% Calculation of Mean and Std
mean_roi = mean(mean(temp_im(temp_im~=0))); % Calculating the mean of the ROI
% THINK ABOUT THE CASE THE IMAGE CONTAINS 0-values!

temp_im = double(temp_im); % convert to double; necessary for the calculation of the std
temp_im_vec = temp_im(:);
std_roi = std(temp_im_vec(temp_im_vec~=0));

lower_limit = 2*std_roi; % lower limit criteria for the edge detection

%% Use criteria to crop ROIs
BW2 = roicolor(image_out(2), mean_roi - lower_limit, max(max(image_out(2))));
labeled_mask = bwlabel(BW2); % labeling of different created masks
ROI_label = labeled_mask.*man_ROI; % only the labeled areas in the inner area
label_max = max(max(labeled_mask)); % Find max label of masks in roi

%% Find ROI
mask_area = 0; % find the label of the mask we really want
if label_max == 0
    warning('There is a problem with the ROI. Please select another ROI.');
end
label = 0; %initialisation
for i = 1:label_max
    mask_area2 = bwarea(ROI_label(ROI_label == i));
    %mask_area2 = bwarea(labeled_mask(labeled_mask == i));
    if mask_area2 > mask_area
        mask_area = mask_area2;
        label = i;
    end
end


%% Fill 'holes' of man_ROI
BW3 = imfill(labeled_mask == label, 'holes');

% %% Cut away parts too far away from the manually selected area
% cropped_area = bwarea(man_ROI); %area of manually selected ROI
% 
labeled_BW = bwlabel(man_ROI);
% 
center = regionprops(labeled_BW,'centroid');
% 
% roi_new = zeros(size(roi,1),2);
% 
% for i = 1:size(roi,1)
% x_dist = roi(i,1)-center.Centroid(1);
% y_dist = roi(i,2)-center.Centroid(2);
% ratio = x_dist/y_dist;
% %% CHANGE THIS LINE IN ORDER TO GET A SMALLER/LARGER CUT-OFF
% new_center_dist = sqrt(x_dist^2+y_dist^2)+ceil((1/1)*sqrt(cropped_area));
% %ORIGINAL CODE begin
% %new_center_dist = sqrt(x_dist^2+y_dist^2)+ceil((1/4)*sqrt(cropped_area));
% %ORIGINAL CODE end
% %%
% if roi(i,2)>center.Centroid(2)
%     y_new = new_center_dist/(sqrt(1+ratio^2));
% else
%     y_new = -new_center_dist/(sqrt(1+ratio^2));
% end
% x_new=y_new*ratio;
% roi_new(i,1)=ceil(roi(i,1)+x_new-x_dist);
% roi_new(i,2)=ceil(roi(i,2)+y_new-y_dist);
% end
% 
% 
% BW_big = roipoly(image_out(2), roi_new(:,1)', roi_new(:,2)'); %create mask

% BW_control = BW_big.*BW3;

BW_control = BW3;

fin_ROI = BW_control;
axis([139,176,267,304])
figure(2)
subplot(1,2,1)
temp_im4 = image_out(2).*double(fin_ROI);
imshow(temp_im4);


%%% BEGINNING NEW ROI CREATION METHOD


figure(3); subplot(1,2,1)
imshow(man_ROI);
subplot(1,2,2)
imshow(fin_ROI);

mean_roi_algorithm = zeros(1,10);
std_roi_algorithm = zeros(1,10);
lower_limit_algorithm = zeros(1,10);
for i=1:10
    temp_im = image_out(2).*BW_control; % crop image to mask area

%     disp('mean:');
    mean_roi_algorithm(i) = mean(mean(temp_im(BW_control~=0))); % Calculating the mean of the ROI

    temp_im = double(temp_im); % convert to double; necessary for the calculation of the std
    temp_im_vec = temp_im(:);
%     disp('std:');
    std_roi_algorithm(i) = std(temp_im_vec(temp_im_vec~=0));

    lower_limit_algorithm(i) = 2*std_roi_algorithm(i); % lower limit criteria for the edge detection
    BW_control = roicolor(image_out(2), mean_roi_algorithm(i) - lower_limit_algorithm(i), max(max(image_out(2))));

    labeled_mask = bwlabel(BW_control); % labeling of different created masks
    ROI_label = labeled_mask.*fin_ROI; % only the labeled areas in the inner area
    label_max = max(max(labeled_mask)); % Find max label of masks in roi


    mask_area = 0; % find the label of the mask we really want
    if label_max == 0
        warning('There is a problem with the ROI. Please select another ROI.');
    end
    label = 0; %initialisation
    for j = 1:label_max
        mask_area2 = bwarea(ROI_label(ROI_label == j));
        %mask_area2 = bwarea(labeled_mask(labeled_mask == i));
        if mask_area2 > mask_area
            mask_area = mask_area2;
            label = j;
        end
    end

    BW3 = imfill(labeled_mask == label, 'holes');


    %BW_control = BW_big.*BW3;
    BW_control = BW3;

    fin_ROI = BW_control;
% 
%     figure
%     temp_im4 = image_out(2).*double(BW_control);
%     imshow(temp_im4);
end

figure(2)
subplot(1,2,2)
temp_im4 = image_out(2).*double(BW_control);
imshow(temp_im4);

%%% END NEW ROI CREATION METHOD

figure(1)
temp_im4 = image_out(2).*double(fin_ROI);
subplot(1,2,2)
imshow(temp_im4);
axis([139,176,267,304])


save 'C:\Documents and Settings\hahntobi\My Documents\MATLAB\temporary saved variables\man_sel_ROI_and_fin_ROI' 'man_ROI' 'mean_roi' 'std_roi' 'BW2' 'BW3' 'center' 'fin_ROI' 'lower_limit' 'mean_roi_algorithm' 'std_roi_algorithm' 'lower_limit_algorithm';