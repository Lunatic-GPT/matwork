function combined_analysis_DCE_study_stand_alone

load det_analysis_gui

FA_ROI_vec = [];
ADC_ROI_vec = [];
ADC_tissue_vec = [];
ADC_tissue2_vec = [];
ADC_lesion_vec = [];
FA_tissue_vec = [];
FA_tissue2_vec = [];
FA_lesion_vec = [];
ADC_tissue_vec_shifted = [];
ADC_tissue2_vec_shifted = [];
ADC_lesion_vec_shifted = [];
FA_tissue_vec_shifted = [];
FA_tissue2_vec_shifted = [];
FA_lesion_vec_shifted = [];
DTI_unw_ROI_vec = [];
DTI_first_ROI_vec = [];
wash_in_ROI_vec = [];
wash_in_tissue_vec = [];
wash_in_tissue2_vec = [];
wash_out_ROI_vec = [];
wash_out_tissue_vec = [];
wash_out_tissue2_vec = [];
diff_lesion_vec = [];
diff_tissue_vec = [];
diff_tissue2_vec = [];
wdiff_lesion_vec = [];
wdiff_tissue_vec = [];
wdiff_tissue2_vec = [];
wash_out_short_ROI_vec = [];
wash_out_short_tissue_vec = [];
wash_out_short_tissue2_vec = [];
DCEimage_tissue_vec = zeros(6,3*512*512);
DCEimage_tissue2_vec = zeros(6,3*512*512);
DCEimage_lesion_vec = zeros(6,3*512*512);
results = cell(1); % result - cell

lesion_volume = 0;

%SID = input('Please enter the patient SID: ','s');
SID = get(analysis,'Enter the subject ID');

    %disp('Please denote now the time differences between the six DCE images [sec]');
    %array_line = input('time point array line: ');
%     array_line = number;
%     load('time_points_array.mat');
%     delta1 = time_points_array(array_line,1);
%     delta2 = time_points_array(array_line,2);
%     delta3 = time_points_array(array_line,3);
%     delta4 = time_points_array(array_line,4);
%     delta5 = time_points_array(array_line,5);

%     shift_x = time_points_array(array_line,6); % THESE ARE THE DCE<->DTI shifts!
%     shift_y = time_points_array(array_line,7);

    shift_x = 0; % we don't need that for this analysis
    shift_y = 0; % "    "

    nslices = get(analysis,'total number of slices');
    first_slice = get(analysis,'first slice');


% nslices = input('How many slices do you want to analyze? (up to 3)');
% 
% 
%                 shift_x = input('Please indicate the shift in the positive x-direction (Matlab convention): ');
%                 shift_y = input('Please indicate the shift in the positive y-direction (Matlab convention): ');

%%% SAVE FOLDER
%save_folder_basic = uigetdir('/export/res/breast','Please select the patient subject folder or the main save folder for the analysis results');
save_folder_basic = ['/export/res/breast/' SID];
%%% END SAVE FOLDER


% cell containing the 'DATA'-folders of the files
load('data_folder.mat');



for slice = 1:nslices
    if slice == 1
%         user_entry_slice = input('Assign the slice number of the DCE series of interest: ');
        user_entry_slice = first_slice;
        user_entry_slice_copy = user_entry_slice;
    else
        user_entry_slice = user_entry_slice+1; % The case that only slice 1 and slice 3 are examined shouldn't happen!
    end
    user_entry_slice_str = num2str(user_entry_slice);
    slicen =num2str(slice);
    save_folder = [save_folder_basic '/research_results_detailed/slice' user_entry_slice_str];
    exist_test = exist(save_folder,'dir');
    if exist_test == 0
        mkdir(save_folder);
    end
    %save_folder2 = [save_folder2_basic '/research_results/slice' user_entry_slice_str];
    %         exist_test = exist(save_folder2,'dir');
    %         if exist_test == 0
    %             mkdir(save_folder2);
    %         end
    close all

    if slice == 1
        %DCE_dir = uigetdir('/home/hahntobi/matlab/Data/','Please select the DCE data folder');
        
%         DCE_dir = data_folder{array_line};

        DCE_dir = get(analysis,'directoryname');
        
        token = strtok(DCE_dir, 'X');
    end
    DCE_dir = [token 'X' user_entry_slice_str];
    final_BW_file = [DCE_dir, '/final_BWs.mat'];
    load(final_BW_file, 'inner_BW_large', 'tissue_BW_large', 'tissue2_BW_large');
    DCE_inner = inner_BW_large;
    DCE_inner_margin = bwmorph(DCE_inner, 'remove');
    DCE_tissue = tissue_BW_large;
    DCE_tissue2 = tissue2_BW_large;
    cut_off_rect_pos = [DCE_dir, '/cut_off_rect_pos.mat'];
    load(cut_off_rect_pos, 'rect_pos');
    x1 = round(rect_pos(1));
    x2 = round(rect_pos(3));
    y1 = round(rect_pos(2));
    y2 = round(rect_pos(4));
    my_axis = [x1 x1+x2 y1 y1+y2];

    DCE_inner_margin_small = DCE_inner_margin((y1:y1+y2-1),x1:(x1+x2-1));


    DCE = zeros(512,512,6); % DCE series array
    DCE_long = zeros(3*x2,6*y2); % all DCE images in one figure
    DCE_relative_change = zeros(x2,6*y2); % seperate image for the relative change images


    source = [DCE_dir, '/ROIs.mat']; %% copy the DCE manually selected area image
    destination = [save_folder, '/manually_selected_ROI_DCE.mat'];
    copyfile(source,destination);
    %destination = [save_folder2, '/manually_selected_ROI_DCE.mat'];
    %copyfile(source,destination);
    copyfile(DCE_dir, save_folder);
    %copyfile(DCE_dir, save_folder2);

    DCE_image_dir = ['/home/hahntobi/matlab/patient_images/' SID '/DCE_temp'];
    for i=1:6 % six DCE images
        DCE_number = 6*(user_entry_slice - 1) + i;
        DCE_number = num2str(DCE_number);
        if (6*(user_entry_slice-1)+i)<10
            DCE_number = ['0' DCE_number];
        end
        file = [DCE_image_dir '/' DCE_number];
        load(file,'image');
        DCE(:,:,i) = image;
    end

    image1 = DCE(:,:,1);
    quot = zeros(6,size(image1,1),size(image1,2)); % relative change images

    %         margin = ones(size(image1,1),size(image1,2));% The gap layer
    %         filled_tissue = imfill(DCE_tissue);
    %         surrounding = bwmorph(filled_tissue, 'remove');
    %         surrounding_small = surrounding((y1:y1+y2-1),x1:(x1+x2-1));
    %         margin = (margin-filled_tissue)+DCE_tissue+DCE_inner;
    %         margin = 1-margin;
    %         margin_small = margin((y1:y1+y2-1),x1:(x1+x2-1));
    %         sur_and_mar = surrounding_small + margin_small;
    %         sur_and_mar_inv = 1-sur_and_mar;

    margin = ones(size(image1,1),size(image1,2));
    filled_tissue2 = imfill(DCE_tissue2);
    surrounding = bwmorph(filled_tissue2, 'remove');
    surrounding_small = surrounding((y1:y1+y2-1),x1:(x1+x2-1));
    margin = (margin-filled_tissue2)+DCE_tissue2+DCE_tissue-bwmorph(imfill(DCE_tissue),'remove')+DCE_inner;
    margin = 1-margin;
    margin_small = margin((y1:y1+y2-1),x1:(x1+x2-1));
    sur_and_mar = surrounding_small + margin_small;
    sur_and_mar_inv = 1-sur_and_mar;

    for i=1:6
        DCEimage = DCE(:,:,i);
        DCE_temp = DCEimage((y1:y1+y2-1),x1:(x1+x2-1));
        DCE_temp_DCE_ROI = DCE_temp.*(1-DCE_inner_margin_small);
        %              DTI_temp_DTI_ROI = DCE_temp.*(1-DTI_inner_margin_small);
        DCE_long(1:x2,((i-1)*x2+1):((i-1)*x2+x2)) = DCE_temp;
        DCE_long(x2+1:2*x2,((i-1)*x2+1):((i-1)*x2+x2)) = DCE_temp_DCE_ROI;
        %              DCE_long(2*x2+1:3*x2,((i-1)*x2+1):((i-1)*x2+x2)) = DTI_temp_DTI_ROI;
        DCE_long(2*x2+1:3*x2,((i-1)*x2+1):((i-1)*x2+x2)) = DCE_temp.*sur_and_mar_inv;
        temp_quot = quot(i,:,:);
        temp_quot = squeeze(temp_quot);
        temp_quot(image1~=0) = (DCEimage(image1~=0)-image1(image1~=0))./image1(image1~=0);
        temp_quot = temp_quot.*(1-DCE_inner_margin); %%!!!!!!!
        DCE_relative_change(1:x2,((i-1)*x2+1):((i-1)*x2+x2)) = temp_quot(y1:(y1+y2-1),x1:(x1+x2-1));
    end

%     figure; imshow(DCE_long, 'InitialMagnification', 'fit'); colormap(jet); colorbar;  caxis([0 2000]);title('DCE series images with DCE tissue and lesion ROI');
% 
%     savename = [save_folder '/DCE_series.fig']; % save as fig...
%     saveas(gcf,savename)
%     savename = [save_folder '/DCE_series.png']; % ... and png.
%     saveas(gcf,savename)

    %         savename = [save_folder2 '/DCE_series.fig']; % save as fig...
    %         saveas(gcf,savename)
    %         savename = [save_folder2 '/DCE_series.png']; % ... and png.
    %         saveas(gcf,savename)
    %
%     figure; imshow(DCE_relative_change, 'InitialMagnification', 'fit'); colormap(jet); colorbar;  caxis([0 2]); title('DCE relative change images');
% 
%     savename = [save_folder '/DCE_rel_change_series.fig']; % save as fig...
%     saveas(gcf,savename)
%     savename = [save_folder '/DCE_rel_change_series.png']; % ... and png.
%     saveas(gcf,savename)

    %         savename = [save_folder2 '/DCE_rel_change_series.fig']; % save as fig...
    %         saveas(gcf,savename)
    %         savename = [save_folder2 '/DCE_rel_change_series.png']; % ... and png.
    %         saveas(gcf,savename)

    %         %% BEGINNING DISPLAY DTI IMAGES & ROI
    %         DTI_filename = ['~hahntobi/matlab/work_images/DTI_temp/slice_' slicen '_avg_DTI'];
    %         load(DTI_filename,'image');
    %         DTI_avg = size_adapt(image);
    %         DTI_avg_small = DTI_avg((y1:y1+y2-1),x1:(x1+x2-1));
    %         fnum = (slice-1)*7+1;
    %         DTI_imagenumber = num2str(fnum);
    %         if fnum < 10
    %             DTI_imagenumber = ['0' DTI_imagenumber];
    %         end
    %         DTI_filename = ['~hahntobi/matlab/work_images/DTI_temp/' DTI_imagenumber];
    %         load(DTI_filename,'image');
    %         unweighted_DTI = size_adapt(image);
    %         unweighted_DTI_small = unweighted_DTI((y1:y1+y2-1),x1:(x1+x2-1));
    %         fnum = (slice-1)*7+2;
    %         DTI_imagenumber = num2str(fnum);
    %         if fnum < 10
    %             DTI_imagenumber = ['0' DTI_imagenumber];
    %         end
    %         DTI_filename = ['~hahntobi/matlab/work_images/DTI_temp/' DTI_imagenumber];
    %         load(DTI_filename,'image');
    %         first_weighted_DTI = size_adapt(image);
    %         first_weighted_DTI_small = first_weighted_DTI((y1:y1+y2-1),x1:(x1+x2-1));
    %         %DTI_long = zeros(x2,4*y2);
    %         DTI_long(:,((1-1)*x2+1):((1-1)*x2+x2)) = unweighted_DTI_small;
    %         DTI_long(:,((2-1)*x2+1):((2-1)*x2+x2)) = first_weighted_DTI_small;
    %         DTI_long(:,((3-1)*x2+1):((3-1)*x2+x2)) = DTI_avg_small;
    %         max_val = max(max(unweighted_DTI_small));
    %         DTI_long(:,((4-1)*x2+1):((4-1)*x2+x2)) = max_val*DTI_inner_margin_small.*ones(x2,y2);
    %         figure; imshow(DTI_long,[]); colorbar
    %         title('DTI images: 1) unweighted 2) 1st weighted 3) avg. DTI 4) ROI from avg. DTI');
    %         savename = [save_folder '/DTI_images.png'];
    %         saveas(gcf,savename)
    %         savename = [save_folder '/DTI_images.fig'];
    %         saveas(gcf,savename)
    %
    %         savename = [save_folder2 '/DTI_images.png'];
    %         saveas(gcf,savename)
    %         savename = [save_folder2 '/DTI_images.fig'];
    %         saveas(gcf,savename)
    %         %% END DISPLAY DTI IMAGES & ROI


    %         % COMPENSATION OF OVERALL SHIFT
    %         DTI_inner_o_shift = shift_image(shift_image(DTI_inner,2,1),2,2); % compensation of overall shift
    %         DTI_inner_b_shift = shift_image(shift_image(DTI_inner,2,1),2,-2); % shifting for best DCE-/DTI-ROI fit
    %         DTI_inner_b_shift_smaller_last_layer = bwmorph(DTI_inner_b_shift, 'remove');
    %         DTI_inner_b_shift_smaller_last_layer_inv = 1-DTI_inner_b_shift_smaller_last_layer;
    %         DTI_inner_b_shift_smaller = DTI_inner_b_shift_smaller_last_layer_inv .* DTI_inner_b_shift;
    %         core_avg = zeros(6,1); % unshifted
    %         rim_avg = zeros(6,1);
    %         core_std = zeros(6,1);
    %         rim_std = zeros(6,1);
    %         %core_o_avg = zeros(6,1);
    %         %rim_o_avg = zeros(6,1);
    %         core_b_avg = zeros(6,1);
    %         rim_b_avg = zeros(6,1);
    %         core_b_std = zeros(6,1);
    %         rim_b_std  = zeros(6,1);
    %         core_b_smaller_avg = zeros(6,1);
    %         rim_b_smaller_avg = zeros(6,1);
    %         core_b_smaller_std = zeros(6,1);
    %         rim_b_smaller_std = zeros(6,1);
    %         tissue_avg = zeros(6,1);
    %         lesion_avg = zeros(6,1);
    %         tissue_std = zeros(6,1);
    %         lesion_std = zeros(6,1);
    DCE_ref = double(DCE(:,:,1));
    DCE_ref(DCE_ref==0)=1; % setting zeros to one, we will forget these pixels later by using the 'no_zero_mask', just set so that the division by DCE_ref is possible
    for i=1:6
        DCE_temp = double(DCE(:,:,i)); % double-conversion to make sure that the calculation of the std is correct
        no_zero_mask = ones(size(DCE_ref,1),size(DCE_ref,2));
        no_zero_mask(DCE_ref==0)=(1-no_zero_mask(DCE_ref==0)); % is zero everywhere where DCE_ref == 0
        no_zero_mask(DCE_temp==0)=0; % If DCE_TEMP is zero at one pixel, this pixel value is wrong!!!!

        %             mask = DTI_inner.*DCE_inner;
        %             temp = no_zero_mask.*((DCE_temp-DCE_ref).*mask)./DCE_ref;
        %             temp_vec = temp(:);
        %             core_avg(i)=mean(mean(temp(mask~=0)));
        %             core_std(i)=std(temp_vec(temp_vec~=0));
        %             mask = (DCE_inner-DTI_inner.*DCE_inner);
        %             temp = no_zero_mask.*((DCE_temp-DCE_ref).*mask)./DCE_ref;
        %             temp_vec = temp(:);
        %             rim_avg(i)=mean(mean(temp(mask~=0)));
        %             rim_std(i)=std(temp_vec(temp_vec~=0));
        %             %mask = DTI_inner_o_shift.*DCE_inner;
        %             %temp = no_zero_mask.*((DCE_temp-DCE_ref).*mask)./DCE_ref;
        %             %core_o_avg(i)=mean(mean(temp(mask~=0)));
        %             %mask = DCE_inner-DTI_inner_o_shift.*DCE_inner;
        %             %temp = no_zero_mask.*((DCE_temp-DCE_ref).*mask)./DCE_ref;
        %             %rim_o_avg(i)=mean(mean(temp(mask~=0)));
        %             mask = DTI_inner_b_shift.*DCE_inner;
        %             temp = no_zero_mask.*((DCE_temp-DCE_ref).*mask)./DCE_ref;
        %             temp_vec = temp(:);
        %             core_b_avg(i)=mean(mean(temp(mask~=0)));
        %             core_b_std(i)=std(temp_vec(temp_vec~=0));
        %             mask = DCE_inner-DTI_inner_b_shift.*DCE_inner;
        %             temp = no_zero_mask.*((DCE_temp-DCE_ref).*mask)./DCE_ref;
        %             temp_vec = temp(:);
        %             rim_b_avg(i)=mean(mean(temp(mask~=0)));
        %             rim_b_std(i)=std(temp_vec(temp_vec~=0));
        %             mask = DTI_inner_b_shift_smaller.*DCE_inner;
        %             temp = no_zero_mask.*((DCE_temp-DCE_ref).*mask)./DCE_ref;
        %             temp_vec = temp(:);
        %             core_b_smaller_avg(i)=mean(mean(temp(mask~=0)));
        %             core_b_smaller_std(i)=std(temp_vec(temp_vec~=0));
        %             mask = DCE_inner-DTI_inner_b_shift_smaller.*DCE_inner;
        %             temp = no_zero_mask.*((DCE_temp-DCE_ref).*mask)./DCE_ref;
        %             temp_vec = temp(:);
        %             rim_b_smaller_avg(i)=mean(mean(temp(mask~=0)));
        %             rim_b_smaller_std(i)=std(temp_vec(temp_vec~=0));


        %% RELATIVE CHANGE CALULATION FOR TIME_SERIES PLOTTING
        mask = DCE_tissue;
        temp = no_zero_mask.*((DCE_temp-DCE_ref).*mask)./DCE_ref;
        temp_vec = temp(mask~=0);
        tissue_avg(i)=mean(temp_vec);
        tissue_std(i)=std(temp_vec);
        mask = DCE_inner;
        temp = no_zero_mask.*((DCE_temp-DCE_ref).*mask)./DCE_ref;
        temp_vec = temp(mask~=0);
        lesion_avg(i)=mean(temp_vec);
        lesion_std(i)=std(temp_vec);
        mask = DCE_tissue2;
        temp = no_zero_mask.*((DCE_temp-DCE_ref).*mask)./DCE_ref;
        temp_vec = temp(mask~=0);
        tissue2_avg(i)=mean(temp_vec);
        tissue2_std(i)=std(temp_vec);
        %% END RELATIVE CHANGE CALCULATION FOR TIME_SERIES

        %% RAW ('BASIC') VALUES CALULATION FOR TIME_SERIES
        mask = DCE_tissue;
        temp = DCE_temp.*mask;
        temp_vec = temp(mask~=0);
        temp_vec(temp_vec==0)=[]; % regard the case that DCE_temp == 0
        r_tissue_avg(i)=mean(temp_vec);
        r_tissue_std(i)=std(temp_vec);
        mask = DCE_inner;
        temp = DCE_temp.*mask;
        temp_vec = temp(mask~=0);
        temp_vec(temp_vec==0)=[];
        r_lesion_avg(i)=mean(temp_vec);
        r_lesion_std(i)=std(temp_vec);
        mask = DCE_tissue2;
        temp = DCE_temp.*mask;
        temp_vec = temp(mask~=0);
        temp_vec(temp_vec==0)=[];
        r_tissue2_avg(i)=mean(temp_vec);
        r_tissue2_std(i)=std(temp_vec);
        %% END RAW VALUES CALCULATION FOR TIME_SERIES

        %% DIFFERENCE CALCULATION FOR TIME_SERIES_DIFF
        mask = DCE_tissue;
        temp = no_zero_mask.*((DCE_temp-DCE_ref).*mask);
        temp_vec = temp(mask~=0);
        d_tissue_avg(i)=mean(temp_vec);
        d_tissue_std(i)=std(temp_vec);
        mask = DCE_inner;
        temp = no_zero_mask.*((DCE_temp-DCE_ref).*mask);
        temp_vec = temp(mask~=0);
        d_lesion_avg(i)=mean(temp_vec);
        d_lesion_std(i)=std(temp_vec);
        mask = DCE_tissue2;
        temp = no_zero_mask.*((DCE_temp-DCE_ref).*mask);
        temp_vec = temp(mask~=0);
        d_tissue2_avg(i)=mean(temp_vec);
        d_tissue2_std(i)=std(temp_vec);
        %% END DIFFERENCE CALCULATION FOR TIME_SERIES_DIFF

        %% WEIGHTED DIFFERENCE CALCULATION FOR TIME_SERIES_WDIFF
        mask = DCE_tissue;
        temp = no_zero_mask.*((DCE_temp-DCE_ref).*mask).*DCE_temp;
        temp_vec = temp(mask~=0);
        wd_tissue_avg(i)=mean(temp_vec);
        wd_tissue_std(i)=std(temp_vec);
        mask = DCE_inner;
        temp = no_zero_mask.*((DCE_temp-DCE_ref).*mask).*DCE_temp;
        temp_vec = temp(mask~=0);
        wd_lesion_avg(i)=mean(temp_vec);
        wd_lesion_std(i)=std(temp_vec);
        mask = DCE_tissue2;
        temp = no_zero_mask.*((DCE_temp-DCE_ref).*mask).*DCE_temp;
        temp_vec = temp(mask~=0);
        wd_tissue2_avg(i)=mean(temp_vec);
        wd_tissue2_std(i)=std(temp_vec);
        %% END WEIGHTED DIFFERENCE CALCULATION FOR TIME_SERIES_WDIFF

        DCEimage_tissue = DCE_temp.*DCE_tissue;
        DCEimage_tissue_vec(i,(user_entry_slice-1)*512*512+1:user_entry_slice*512*512) = DCEimage_tissue(:)';
        if slice == nslices
            temp = DCEimage_tissue_vec(i,:);
            DCEimage_tissue_mean(i) = mean(temp(temp~=0));
            DCEimage_tissue_std(i) = std(temp(temp~=0));
        end

        DCEimage_tissue2 = DCE_temp.*DCE_tissue2;
        DCEimage_tissue2_vec(i,(user_entry_slice-1)*512*512+1:user_entry_slice*512*512) = DCEimage_tissue2(:)';
        if slice == nslices
            temp = DCEimage_tissue2_vec(i,:);
            DCEimage_tissue2_mean(i) = mean(temp(temp~=0));
            DCEimage_tissue2_std(i) = std(temp(temp~=0));
        end

        DCEimage_lesion = DCE_temp.*DCE_inner;
        DCEimage_lesion_vec(i,(user_entry_slice-1)*512*512+1:user_entry_slice*512*512) = DCEimage_lesion(:)';
        if slice == nslices
            temp = DCEimage_lesion_vec(i,:);
            DCEimage_lesion_mean(i) = mean(temp(temp~=0));
            DCEimage_lesion_std(i) = std(temp(temp~=0));
        end
        %         figure; subplot(3,2,1); imshow(DCE(:,:,i).*DTI_inner_o_shift.*DCE_inner,[]);
        %         subplot(3,2,2); imshow(DCE(:,:,i).*(DCE_inner-DTI_inner_o_shift.*DCE_inner),[]);
        %         subplot(3,2,3); imshow(DCE(:,:,i).*DTI_inner_b_shift.*DCE_inner,[]);
        %         subplot(3,2,4); imshow(DCE(:,:,i).*(DCE_inner-DTI_inner_b_shift.*DCE_inner),[]);
        %         subplot(3,2,5); imshow(DCE(:,:,i).*DCE_tissue,[]);
    end

    %         results = [results core_b_avg 'core_b_avg' core_std 'core_std' rim_b_avg 'rim_b_avg' rim_std 'rim_std' core_b_smaller_avg 'core_b_smaller_avg' core_b_smaller_std 'core_b_smaller_std' rim_b_smaller_avg 'rim_b_smaller_avg' rim_b_smaller_std 'rim_b_smaller_std' tissue_avg 'tissue_avg' tissue_std 'tissue_std' lesion_avg 'lesion_avg' lesion_std 'lesion_std'];

    %%% ADC image analysis start
%     file = ['/home/hahntobi/matlab/patient_images/' SID '/ADC_temp/0' user_entry_slice_str '.mat'];
%     load(file,'image');
%     ADCimage = image;
%     ADCimage = size_adapt(ADCimage);
% 
%     ADC_temp = ADCimage((y1:y1+y2-1),x1:(x1+x2-1));
%     %ADC_temp_DCE_ROI = ADC_temp.*(1-DCE_inner_margin_small); % should
%     %be the same as below
%     ADC_temp_DCE_ROI = ADC_temp.*sur_and_mar_inv;
% %     figure; imshow(ADC_temp_DCE_ROI, 'InitialMagnification', 'fit'); colorbar; caxis([0 3]); %colormap(jet)
% % 
% %     savename = [save_folder '/ADC_with_DCE_ROI.fig']; % save as fig...
% %     saveas(gcf,savename)
% %     savename = [save_folder '/ADC_with_DCE_ROI.png']; % ... and png.
% %     saveas(gcf,savename)
% 
%     %         savename = [save_folder2 '/ADC_with_DCE_ROI.fig']; % save as fig...
%     %         saveas(gcf,savename)
%     %         savename = [save_folder2 '/ADC_with_DCE_ROI.png']; % ... and png.
%     %         saveas(gcf,savename)
% 
%     %% ADC analysis with DCE ROI
%     ADC_ROI = ADCimage.*DCE_inner;
%     ADC_ROI_vec = [ADC_ROI_vec ADC_ROI(:)];
%     if slice == nslices % plot histogram after last slice
% %         figure; hist(ADC_ROI_vec(ADC_ROI_vec~=0),50);
% %         savename = [save_folder '/ADC_hist.png'];
% %         saveas(gcf,savename)
% %         savename = [save_folder '/ADC_hist.fig'];
% %         saveas(gcf,savename)
% 
%         %             savename = [save_folder2 '/ADC_hist.png'];
%         %             saveas(gcf,savename)
%         %             savename = [save_folder2 '/ADC_hist.fig'];
%         %             saveas(gcf,savename)
%         %ADC_ROI_DCE_mean = mean(ADC_ROI_vec(ADC_ROI_vec~=0));
%         %ADC_ROI_DCE_std = std(ADC_ROI_vec(ADC_ROI_vec~=0));
%     end
% 
%     ADC_tissue = ADCimage.*DCE_tissue;
%     ADC_tissue_vec = [ADC_tissue_vec ADC_tissue(:)];
%     if slice == nslices
%         ADC_tissue_mean = mean(ADC_tissue_vec(ADC_tissue_vec~=0));
%         ADC_tissue_std = std(ADC_tissue_vec(ADC_tissue_vec~=0));
%     end
% 
%     ADC_tissue2 = ADCimage.*DCE_tissue2;
%     ADC_tissue2_vec = [ADC_tissue2_vec ADC_tissue2(:)];
%     if slice == nslices
%         ADC_tissue2_mean = mean(ADC_tissue2_vec(ADC_tissue2_vec~=0));
%         ADC_tissue2_std = std(ADC_tissue2_vec(ADC_tissue2_vec~=0));
%     end
% 
%     ADC_lesion = ADCimage.*DCE_inner;
%     ADC_lesion_vec = [ADC_lesion_vec ADC_lesion(:)];
%     if slice == nslices
%         ADC_lesion_mean = mean(ADC_lesion_vec(ADC_lesion_vec~=0));
%         ADC_lesion_std = std(ADC_lesion_vec(ADC_lesion_vec~=0));
%     end
% 
%     %% END ADC analysis with DCE ROI
% 
%     %             %% ADC analysis with DTI ROI
%     %             ADC_ROI = ADCimage.*DTI_inner;
%     %             ADC_ROI_vec = [ADC_ROI_vec ADC_ROI(:)];
%     %             if slice == 3 % plot histogram after third slice
%     %                 figure; hist(ADC_ROI_vec(ADC_ROI_vec~=0),50);
%     %                 savename = [save_folder '/ADC_hist.png'];
%     %                 saveas(gcf,savename)
%     %                 savename = [save_folder '/ADC_hist.fig'];
%     %                 saveas(gcf,savename)
%     %
%     %                 savename = [save_folder2 '/ADC_hist.png'];
%     %                 saveas(gcf,savename)
%     %                 savename = [save_folder2 '/ADC_hist.fig'];
%     %                 saveas(gcf,savename)
%     %             end
%     %             figure; subplot(2,2,1);imshow(ADC_ROI,[]);
%     %             axis(my_axis); title('DTI ROI assigned to the ADC image');
%     %             my_caxis = caxis;
%     %             temp = ADC_ROI(:);
%     %             temp = temp(temp~=0);
%     %             ADC_ROI_core=ADC_ROI;
%     %             ADC_ROI_core(ADC_ROI_core<mean(temp))=0;
%     %             ADC_ROI_rim = ADC_ROI - ADC_ROI_core;
%     %             ADC_mean_core = mean(mean(ADC_ROI_core(ADC_ROI_core~=0)));
%     %             ADC_mean_rim = mean(mean(ADC_ROI_rim(ADC_ROI_rim~=0)));
%     %             subplot(2,2,2);imshow(ADC_ROI_core,my_caxis);axis(my_axis); title('Core and outer part ADC image');
%     %             subplot(2,2,3);imshow(ADC_ROI_rim, my_caxis); axis(my_axis); title('Rim ADC image');
%     %             savename = [save_folder '/ADC_with_ROI.png'];
%     %             saveas(gcf,savename)
%     %             savename = [save_folder '/ADC_with_ROI.fig'];
%     %             saveas(gcf,savename)
%     %
%     %             savename = [save_folder2 '/ADC_with_ROI.png'];
%     %             saveas(gcf,savename)
%     %             savename = [save_folder2 '/ADC_with_ROI.fig'];
%     %             saveas(gcf,savename)
%     %             %% END ADC analysis with DTI ROI
% 
%     %%% ADC image analysis end

    %%% Difference % Weighted Difference image analysis start
    file = ['/home/hahntobi/matlab/patient_images/' SID '/DCE_temp/diff_slice' user_entry_slice_str '.mat'];
    load(file,'image');
    diffimage = image;
    diffimage_temp = diffimage((y1:y1+y2-1),x1:(x1+x2-1));
    %diffimage_DCE_ROI = diffimage_temp.*sur_and_mar_inv; % WITH
    %MARGINS
    diffimage_DCE_ROI = diffimage_temp;
    figure; imshow(diffimage_DCE_ROI, 'InitialMagnification', 'fit'); colorbar; colormap(jet);caxis([0 1600]);%colormap(jet);
    savename = [save_folder '/diff.fig']; % save as fig...
    saveas(gcf,savename)
    savename = [save_folder '/diff.png']; % ... and png.
    saveas(gcf,savename)

    wdiffimage = diffimage.*DCE(:,:,2); % weighted Difference image
    diffweight_temp = wdiffimage((y1:y1+y2-1),x1:(x1+x2-1));
    %diffweight_DCE_ROI = diffweight_temp.*sur_and_mar_inv; % WITH
    %MARGINS
    diffweight_DCE_ROI = diffweight_temp;
%     figure; imshow(diffweight_DCE_ROI, 'InitialMagnification', 'fit'); colorbar; colormap(jet); caxis([0 2E6]);%colormap(jet);
%     savename = [save_folder '/wdiff.fig']; % save as fig...
%     saveas(gcf,savename)
%     savename = [save_folder '/wdiff.png']; % ... and png.
%     saveas(gcf,savename)

    diff_tissue = diffimage.*DCE_tissue;
    diff_tissue_vec = [diff_tissue_vec diff_tissue(:)];
    if slice == nslices
        diff_tissue_mean = mean(diff_tissue_vec(diff_tissue_vec~=0));
        diff_tissue_std = std(diff_tissue_vec(diff_tissue_vec~=0));
    end

    diff_tissue2 = diffimage.*DCE_tissue2;
    diff_tissue2_vec = [diff_tissue2_vec diff_tissue2(:)];
    if slice == nslices
        diff_tissue2_mean = mean(diff_tissue2_vec(diff_tissue2_vec~=0));
        diff_tissue2_std = std(diff_tissue2_vec(diff_tissue2_vec~=0));
    end

    diff_lesion = diffimage.*DCE_inner;
    diff_lesion_vec = [diff_lesion_vec diff_lesion(:)];
    if slice == nslices
        diff_lesion_mean = mean(diff_lesion_vec(diff_lesion_vec~=0));
        diff_lesion_std = std(diff_lesion_vec(diff_lesion_vec~=0));
    end

    wdiff_tissue = wdiffimage.*DCE_tissue;
    wdiff_tissue_vec = [wdiff_tissue_vec wdiff_tissue(:)];
    if slice == nslices
        wdiff_tissue_mean = mean(wdiff_tissue_vec(wdiff_tissue_vec~=0));
        wdiff_tissue_std = std(wdiff_tissue_vec(wdiff_tissue_vec~=0));
    end

    wdiff_tissue2 = wdiffimage.*DCE_tissue2;
    wdiff_tissue2_vec = [wdiff_tissue2_vec wdiff_tissue2(:)];
    if slice == nslices
        wdiff_tissue2_mean = mean(wdiff_tissue2_vec(wdiff_tissue2_vec~=0));
        wdiff_tissue2_std = std(wdiff_tissue2_vec(wdiff_tissue2_vec~=0));
    end

    wdiff_lesion = wdiffimage.*DCE_inner;
    wdiff_lesion_vec = [wdiff_lesion_vec wdiff_lesion(:)];
    if slice == nslices
        wdiff_lesion_mean = mean(wdiff_lesion_vec(wdiff_lesion_vec~=0));
        wdiff_lesion_std = std(wdiff_lesion_vec(wdiff_lesion_vec~=0));
    end
    %%% Difference % Weighted Difference image analysis end


    %%% FA image analysis start
    %         [filename,pathname] = uigetfile({'*.m;*.mat','MATLAB Files (*.m,*.mat)'},'Select the slice corresponding FA image','/home/hahntobi/matlab/patient_images/');
    %         file = [pathname filename];
%     file = ['/home/hahntobi/matlab/patient_images/' SID '/FA_temp/0' user_entry_slice_str '.mat'];
%     load(file,'image');
%     FAimage = image;
%     FAimage = size_adapt(FAimage);
% 
%     FA_temp = FAimage((y1:y1+y2-1),x1:(x1+x2-1));
%     %FA_temp_DCE_ROI = FA_temp.*(1-DCE_inner_margin_small);
%     FA_temp_DCE_ROI = FA_temp.*sur_and_mar_inv;
% %     figure; imshow(FA_temp_DCE_ROI, 'InitialMagnification', 'fit'); colorbar; caxis([0 1.2]);%colormap(jet);
% % 
% %     savename = [save_folder '/FA_with_DCE_ROI.fig']; % save as fig...
% %     saveas(gcf,savename)
% %     savename = [save_folder '/FA_with_DCE_ROI.png']; % ... and png.
% %     saveas(gcf,savename)
%     %
%     %         savename = [save_folder2 '/FA_with_DCE_ROI.fig']; % save as fig...
%     %         saveas(gcf,savename)
%     %         savename = [save_folder2 '/FA_with_DCE_ROI.png']; % ... and png.
%     %         saveas(gcf,savename)
% 
%     %% FA analysis with DCE ROI
%     FA_ROI = FAimage.*DCE_inner;
%     FA_ROI_vec = [FA_ROI_vec FA_ROI(:)];
%     if slice == nslices % plot histogram after last slice
% %         figure; hist(FA_ROI_vec(FA_ROI_vec~=0),50);
% %         savename = [save_folder '/FA_hist.png'];
% %         saveas(gcf,savename)
% %         savename = [save_folder '/FA_hist.fig'];
% %         saveas(gcf,savename)
% 
%         %             savename = [save_folder2 '/FA_hist.png'];
%         %             saveas(gcf,savename)
%         %             savename = [save_folder2 '/FA_hist.fig'];
%         %             saveas(gcf,savename)
%         %FA_ROI_DCE_mean = mean(FA_ROI_vec(FA_ROI_vec~=0));
%         %FA_ROI_DCE_std = std(FA_ROI_vec(FA_ROI_vec~=0));
%     end
% 
%     FA_tissue = FAimage.*DCE_tissue;
%     FA_tissue_vec = [FA_tissue_vec FA_tissue(:)];
%     if slice == nslices
%         FA_tissue_mean = mean(FA_tissue_vec(FA_tissue_vec~=0));
%         FA_tissue_std = std(FA_tissue_vec(FA_tissue_vec~=0));
%     end
% 
%     FA_tissue2 = FAimage.*DCE_tissue2;
%     FA_tissue2_vec = [FA_tissue2_vec FA_tissue2(:)];
%     if slice == nslices
%         FA_tissue2_mean = mean(FA_tissue2_vec(FA_tissue2_vec~=0));
%         FA_tissue2_std = std(FA_tissue2_vec(FA_tissue2_vec~=0));
%     end
% 
%     FA_lesion = FAimage.*DCE_inner;
%     FA_lesion_vec = [FA_lesion_vec FA_lesion(:)];
%     if slice == nslices
%         FA_lesion_mean = mean(FA_lesion_vec(FA_lesion_vec~=0));
%         FA_lesion_std = std(FA_lesion_vec(FA_lesion_vec~=0));
%     end
%     %% END FA analysis with DCE ROI
% 
%     %% DTI analysis with DCE ROI
%     %         [filename,pathname] = uigetfile({'*.m;*.mat','MATLAB Files (*.m,*.mat)'},'Select the slice corresponding unweighted DTI image','/home/hahntobi/matlab/patient_images/');
%     %         file = [pathname filename];
%     DTI_number = (user_entry_slice-1)*7+1;
%     DTI_number_str = num2str(DTI_number);
%     if user_entry_slice > 2
%         file = ['/home/hahntobi/matlab/patient_images/' SID '/DTI_temp/' DTI_number_str '.mat'];
%     else
%         file = ['/home/hahntobi/matlab/patient_images/' SID '/DTI_temp/0' DTI_number_str '.mat'];
%     end
%     load(file,'image');
%     unw_DTI = image;
%     unw_DTI = size_adapt(unw_DTI);
%     unw_DTI_temp = unw_DTI((y1:y1+y2-1),x1:(x1+x2-1));
%     unw_DTI_temp_DCE_ROI = unw_DTI_temp.*(1-DCE_inner_margin_small);
% %     figure; imshow(unw_DTI_temp_DCE_ROI, 'InitialMagnification', 'fit'); colorbar; caxis([0 800]); %colormap(jet);
% % 
% %     savename = [save_folder '/unw_DTI_with_DCE_ROI.fig']; % save as fig...
% %     saveas(gcf,savename)
% %     savename = [save_folder '/unw_DTI_with_DCE_ROI.png']; % ... and png.
% %     saveas(gcf,savename)
% 
%     %         savename = [save_folder2 '/unw_DTI_with_DCE_ROI.fig']; % save as fig...
%     %         saveas(gcf,savename)
%     %         savename = [save_folder2 '/unw_DTI_with_DCE_ROI.png']; % ... and png.
%     %         saveas(gcf,savename)
%     %
%     %         [filename,pathname] = uigetfile({'*.m;*.mat','MATLAB Files (*.m,*.mat)'},'Select the slice corresponding first weighted DTI image','/home/hahntobi/matlab/patient_images/');
%     %         file = [pathname filename];
%     DTI_number = (user_entry_slice-1)*7+2;
%     DTI_number_str = num2str(DTI_number);
%     if user_entry_slice > 2
%         file = ['/home/hahntobi/matlab/patient_images/' SID '/DTI_temp/' DTI_number_str '.mat'];
%     else
%         file = ['/home/hahntobi/matlab/patient_images/' SID '/DTI_temp/0' DTI_number_str '.mat'];
%     end
%     load(file,'image');
%     first_DTI = image;
%     first_DTI = size_adapt(first_DTI);
% 
%     first_DTI_temp = first_DTI((y1:y1+y2-1),x1:(x1+x2-1));
%     first_DTI_temp_DCE_ROI = first_DTI_temp.*(1-DCE_inner_margin_small);
% %     figure; imshow(first_DTI_temp_DCE_ROI, 'InitialMagnification', 'fit'); colorbar; caxis([0 400]);%colormap(jet);
% % 
% %     savename = [save_folder '/first_DTI_with_DCE_ROI.fig']; % save as fig...
% %     saveas(gcf,savename)
% %     savename = [save_folder '/first_DTI_with_DCE_ROI.png']; % ... and png.
% %     saveas(gcf,savename)
%     %
%     %         savename = [save_folder2 '/first_DTI_with_DCE_ROI.fig']; % save as fig...
%     %         saveas(gcf,savename)
%     %         savename = [save_folder2 '/first_DTI_with_DCE_ROI.png']; % ... and png.
%     %         saveas(gcf,savename)
% 
%     DTI_unw_ROI = unw_DTI.*DCE_inner;
%     DTI_unw_ROI_vec = [DTI_unw_ROI_vec DTI_unw_ROI(DCE_inner~=0)'];
%     if slice == nslices
%         DTI_unw_ROI_mean = mean(DTI_unw_ROI_vec);
%         DTI_unw_ROI_std = std(DTI_unw_ROI_vec);
%     end
% 
%     DTI_first_ROI = first_DTI.*DCE_inner;
%     DTI_first_ROI_vec = [DTI_first_ROI_vec DTI_first_ROI(:)];
%     if slice == nslices
%         DTI_first_ROI_mean = mean(DTI_first_ROI_vec(DTI_first_ROI_vec~=0));
%         DTI_first_ROI_std = std(DTI_first_ROI_vec(DTI_first_ROI_vec~=0));
%     end
% 
% 
%     file = ['/home/hahntobi/matlab/patient_images/' SID '/DTI_temp/' 'slice_' user_entry_slice_str '_avg_DTI.mat'];
%     load(file,'image');
%     avg_DTI = image;
%     avg_DTI = size_adapt(avg_DTI);
% 
%     avg_DTI_temp = avg_DTI((y1:y1+y2-1),x1:(x1+x2-1));
%     avg_DTI_temp_DCE_ROI = avg_DTI_temp.*(1-DCE_inner_margin_small);
% %     figure; imshow(avg_DTI_temp_DCE_ROI, 'InitialMagnification', 'fit'); colorbar; caxis([0 400]);%colormap(jet);
% % 
% %     savename = [save_folder '/avg_DTI_with_DCE_ROI.fig']; % save as fig...
% %     saveas(gcf,savename)
% %     savename = [save_folder '/avg_DTI_with_DCE_ROI.png']; % ... and png.
% %     saveas(gcf,savename)
% 
%     %         savename = [save_folder2 '/avg_DTI_with_DCE_ROI.fig']; % save as fig...
%     %         saveas(gcf,savename)
%     %         savename = [save_folder2 '/avg_DTI_with_DCE_ROI.png']; % ... and png.
%     %         saveas(gcf,savename)

    %% END DTI analysis with DCE ROI
    %
    %             %% FA analysis with DTI ROI
    %             FA_ROI = FAimage.*DTI_inner;
    %             FA_ROI_vec = [FA_ROI_vec FA_ROI(:)];
    %             if slice == 3 % plot histogram after third slice
    %                 figure; hist(FA_ROI_vec(FA_ROI_vec~=0),50);
    %                 savename = [save_folder '/FA_hist.png'];
    %                 saveas(gcf,savename)
    %                 savename = [save_folder '/FA_hist.fig'];
    %                 saveas(gcf,savename)
    %
    %                 savename = [save_folder2 '/FA_hist.png'];
    %                 saveas(gcf,savename)
    %                 savename = [save_folder2 '/FA_hist.fig'];
    %                 saveas(gcf,savename)
    %             end
    %             figure; subplot(2,2,1);imshow(FA_ROI,[]);
    %             axis(my_axis); title('DTI ROI assigned to the FA image');
    %             my_caxis = caxis;
    %             temp = FA_ROI(:);
    %             temp = temp(temp~=0);
    %             FA_ROI_core=FA_ROI;
    %             FA_ROI_core(FA_ROI_core<mean(temp))=0;
    %             FA_ROI_rim = FA_ROI - FA_ROI_core;
    %             FA_mean_core = mean(mean(FA_ROI_core(FA_ROI_core~=0)));
    %             FA_mean_rim = mean(mean(FA_ROI_rim(FA_ROI_rim~=0)));
    %             subplot(2,2,2);imshow(FA_ROI_core,my_caxis);axis(my_axis); title('Core and outer part FA image');
    %             subplot(2,2,3);imshow(FA_ROI_rim, my_caxis); axis(my_axis); title('Rim FA image');
    %             savename = [save_folder '/FA_with_ROI.png'];
    %             saveas(gcf,savename)
    %             savename = [save_folder '/FA_with_ROI.fig'];
    %             saveas(gcf,savename)
    %
    %             savename = [save_folder2 '/FA_with_ROI.png'];
    %             saveas(gcf,savename)
    %             savename = [save_folder2 '/FA_with_ROI.fig'];
    %             saveas(gcf,savename)
    %             %% END FA analysis with DTI ROI
    %%% FA image analysis end


    %% WASH-IN Analysis
    %         [filename,pathname] = uigetfile({'*.m;*.mat','MATLAB Files (*.m,*.mat)'},'Select the slice corresponding WASH-IN image','/home/hahntobi/matlab/patient_images/');
    %         file = [pathname filename];
    file = ['/home/hahntobi/matlab/patient_images/' SID '/DCE_temp/wash_in_slice' user_entry_slice_str '.mat'];
    load(file,'image');
    wash_in = image;
    wash_in_ROI = wash_in.*DCE_inner;
    wash_in_ROI_single_vec = wash_in_ROI(:);
%     figure; hist(wash_in_ROI_single_vec(wash_in_ROI_single_vec~=0),50);
%     savename = [save_folder '/wash_in_hist.png'];
%     saveas(gcf,savename)
%     savename = [save_folder '/wash_in_hist.fig'];
%     saveas(gcf,savename)
    %         savename = [save_folder2 '/wash_in_hist.png'];
    %         saveas(gcf,savename)
    %         savename = [save_folder2 '/wash_in_hist.fig'];
    %         saveas(gcf,savename)


    wash_in_ROI_vec = [wash_in_ROI_vec wash_in_ROI_single_vec];
    if slice == nslices
        wash_in_ROI_DCE_mean = mean(wash_in_ROI_vec(wash_in_ROI_vec~=0));
        wash_in_ROI_DCE_std = std(wash_in_ROI_vec(wash_in_ROI_vec~=0));
    end

    wash_in_tissue = wash_in.*DCE_tissue;
    wash_in_tissue_vec = [wash_in_tissue_vec wash_in_tissue(:)];
    if slice == nslices
        wash_in_tissue_mean = mean(wash_in_tissue_vec(wash_in_tissue_vec~=0));
        wash_in_tissue_std = std(wash_in_tissue_vec(wash_in_tissue_vec~=0));
    end

    wash_in_tissue2 = wash_in.*DCE_tissue2;
    wash_in_tissue2_vec = [wash_in_tissue2_vec wash_in_tissue2(:)];
    if slice == nslices
        wash_in_tissue2_mean = mean(wash_in_tissue2_vec(wash_in_tissue2_vec~=0));
        wash_in_tissue2_std = std(wash_in_tissue2_vec(wash_in_tissue2_vec~=0));
    end

    wash_in_temp  = wash_in((y1:y1+y2-1),x1:(x1+x2-1)); % = DCE_long(3*x2+1:4*x2,((2-1)*x2+1):((2-1)*x2+x2)) SHOULD BE THE SAME!!!
%     figure; imshow(wash_in_temp, 'InitialMagnification', 'fit'); colormap(jet);colorbar; caxis([0 2])
% 
%     savename = [save_folder '/wash_in.fig']; % save as fig...
%     saveas(gcf,savename)
%     savename = [save_folder '/wash_in.png']; % ... and png.
%     saveas(gcf,savename)

    %         savename = [save_folder2 '/wash_in.fig']; % save as fig...
    %         saveas(gcf,savename)
    %         savename = [save_folder2 '/wash_in.png']; % ... and png.
    %         saveas(gcf,savename)
    %% END WASH-IN Analysis

    %% WASH-OUT Analysis
    file = ['/home/hahntobi/matlab/patient_images/' SID '/DCE_temp/wash_out' user_entry_slice_str '.mat'];
    load(file,'image');
    wash_out = image;
    wash_out_ROI = wash_out.*DCE_inner;
    wash_out_ROI_single_vec = wash_out_ROI(:);
%     figure; hist(wash_out_ROI_single_vec(wash_out_ROI_single_vec~=0),50);
%     savename = [save_folder '/wash_out_hist.png'];
%     saveas(gcf,savename)
%     savename = [save_folder '/wash_out_hist.fig'];
%     saveas(gcf,savename)
    %         savename = [save_folder2 '/wash_out_hist.png'];
    %         saveas(gcf,savename)
    %         savename = [save_folder2 '/wash_out_hist.fig'];
    %         saveas(gcf,savename)

    wash_out_ROI_vec = [wash_out_ROI_vec wash_out_ROI_single_vec];
    if slice == nslices
        wash_out_ROI_DCE_mean = mean(wash_out_ROI_vec(wash_out_ROI_vec~=0));
        wash_out_ROI_DCE_std = std(wash_out_ROI_vec(wash_out_ROI_vec~=0));
    end

    wash_out_tissue = wash_out.*DCE_tissue;
    wash_out_tissue_vec = [wash_out_tissue_vec wash_out_tissue(:)];
    if slice == nslices
        wash_out_tissue_mean = mean(wash_out_tissue_vec(wash_out_tissue_vec~=0));
        wash_out_tissue_std = std(wash_out_tissue_vec(wash_out_tissue_vec~=0));
    end

    wash_out_tissue2 = wash_out.*DCE_tissue2;
    wash_out_tissue2_vec = [wash_out_tissue2_vec wash_out_tissue2(:)];
    if slice == nslices
        wash_out_tissue2_mean = mean(wash_out_tissue2_vec(wash_out_tissue2_vec~=0));
        wash_out_tissue2_std = std(wash_out_tissue2_vec(wash_out_tissue2_vec~=0));
    end

    wash_out_temp  = wash_out((y1:y1+y2-1),x1:(x1+x2-1)); % = DCE_long(3*x2+1:4*x2,((2-1)*x2+1):((2-1)*x2+x2)) SHOULD BE THE SAME!!!
%     figure; imshow(wash_out_temp, 'InitialMagnification', 'fit')
%     caxis([-90 90])
%     jet_inv = jet;
%     jet_inv = flipud(jet_inv);
%     colormap(jet_inv)
%     colorbar
% 
%     savename = [save_folder '/wash_out.fig']; % save as fig...
%     saveas(gcf,savename)
%     savename = [save_folder '/wash_out.png']; % ... and png.
%     saveas(gcf,savename)

    %         savename = [save_folder2 '/wash_out.fig']; % save as fig...
    %         saveas(gcf,savename)
    %         savename = [save_folder2 '/wash_out.png']; % ... and png.
    %         saveas(gcf,savename)
    %% END WASH-OUT Analysis

    %% 'SHORT' WASH-OUT Analysis
    file = ['/home/hahntobi/matlab/patient_images/' SID '/DCE_temp/wash_out_short' user_entry_slice_str '.mat'];
    load(file,'image');
    wash_out_short = image;
    wash_out_short_ROI = wash_out_short.*DCE_inner;
    wash_out_short_ROI_single_vec = wash_out_short_ROI(:);
%     figure; hist(wash_out_short_ROI_single_vec(wash_out_short_ROI_single_vec~=0),50);
%     savename = [save_folder '/wash_out_short_hist.png'];
%     saveas(gcf,savename)
%     savename = [save_folder '/wash_out_short_hist.fig'];
%     saveas(gcf,savename)
    %         savename = [save_folder2 '/wash_out_short_hist.png'];
    %         saveas(gcf,savename)
    %         savename = [save_folder2 '/wash_out_short_hist.fig'];
    %         saveas(gcf,savename)

    wash_out_short_ROI_vec = [wash_out_short_ROI_vec wash_out_short_ROI_single_vec];
    if slice == nslices
        wash_out_short_ROI_DCE_mean = mean(wash_out_short_ROI_vec(wash_out_short_ROI_vec~=0));
        wash_out_short_ROI_DCE_std = std(wash_out_short_ROI_vec(wash_out_short_ROI_vec~=0));
    end

    wash_out_short_tissue = wash_out_short.*DCE_tissue;
    wash_out_short_tissue_vec = [wash_out_short_tissue_vec wash_out_short_tissue(:)];
    if slice == nslices
        wash_out_short_tissue_mean = mean(wash_out_short_tissue_vec(wash_out_short_tissue_vec~=0));
        wash_out_short_tissue_std = std(wash_out_short_tissue_vec(wash_out_short_tissue_vec~=0));
    end

    wash_out_short_tissue2 = wash_out_short.*DCE_tissue2;
    wash_out_short_tissue2_vec = [wash_out_short_tissue2_vec wash_out_short_tissue2(:)];
    if slice == nslices
        wash_out_short_tissue2_mean = mean(wash_out_short_tissue2_vec(wash_out_short_tissue2_vec~=0));
        wash_out_short_tissue2_std = std(wash_out_short_tissue2_vec(wash_out_short_tissue2_vec~=0));
    end

    wash_out_short_temp  = wash_out_short((y1:y1+y2-1),x1:(x1+x2-1)); % = DCE_long(3*x2+1:4*x2,((2-1)*x2+1):((2-1)*x2+x2)) SHOULD BE THE SAME!!!
%     figure; imshow(wash_out_short_temp, 'InitialMagnification', 'fit')
%     caxis([-90 90])
%     jet_inv = jet;
%     jet_inv = flipud(jet_inv);
%     colormap(jet_inv)
%     colorbar
% 
%     savename = [save_folder '/wash_out_short.fig']; % save as fig...
%     saveas(gcf,savename)
%     savename = [save_folder '/wash_out_short.png']; % ... and png.
%     saveas(gcf,savename)

    %         savename = [save_folder2 '/wash_out_short.fig']; % save as fig...
    %         saveas(gcf,savename)
    %         savename = [save_folder2 '/wash_out_short.png']; % ... and png.
    %         saveas(gcf,savename)
    %% END 'SHORT' WASH-OUT Analysis


    disp('____________________________________________________________________');
    disp('                      ANALYSIS RESULTS                              ');
    %         pixels_outside_o = sum(sum((DTI_inner_o_shift-DCE_inner)==1));
    %         pixels_outside_o = num2str(pixels_outside_o);
    %         message = ['There are ', pixels_outside_o, ' pixels of the overall shift compensated DTI ROI outside of the DCE ROI.'];
    %         save_text = [save_text message;];
    %         disp(message);
    %         pixels_outside_b = sum(sum((DTI_inner_b_shift-DCE_inner)==1));
    %         pixels_outside_b_str = num2str(pixels_outside_b);
    %         results = [results pixels_outside_b 'pixels of the optimal shifted DTI ROI outside of the DCE ROI.'];
    %         message = ['There are ', pixels_outside_b_str, ' pixels of the optimal shifted DTI ROI outside of the DCE ROI.'];
    %         save_text = [save_text message];
    %         disp(message);
    DCE_ROI_area = [DCE_dir, '/ROIs.mat'];
    load(DCE_ROI_area,'fin_area');
    lesion_volume = lesion_volume + fin_area;
    if slice == nslices
        results = [results lesion_volume 'lesion volume']; % we area really only interested in the volume of the whole lesion, not in the area of the single slices
    end
%         results = [results fin_area 'area of the DCE ROI'];
    %         DTI_ROI_area = [DTI_dir, '/ROIs.mat'];
    %         load(DTI_ROI_area,'fin_area');
    %         fin_area_str = num2str(4*fin_area);
    %         message = ['The area of the DTI ROI is ' fin_area_str ' pixels.'];
    %         save_text = [save_text message];
    %         results = [results fin_area 'area of the DTI ROI'];
    %         disp(message);
    %         load(DCE_ROI_area,'lower_limit', 'std_roi');
    %         lower_limit_roi = lower_limit/std_roi;
    %         lower_limit_str = num2str(lower_limit);
    %         lower_limit_roi_str = num2str(lower_limit_roi);
    %         message = ['The lower limit criterion for the DCE ROI: ' lower_limit_str ', i.e. ' lower_limit_roi_str ' * std(manually selected area)'];
    %         save_text = [save_text message];
    %         results = [results lower_limit 'The lower limit criterion for the DCE ROI' lower_limit_roi '*std(manually selected area)'];
    %         disp(message);
    %         load(DTI_ROI_area,'lower_limit', 'std_roi');
    %         lower_limit_roi = lower_limit/std_roi;
    %         lower_limit_str = num2str(lower_limit);
    %         lower_limit_roi_str = num2str(lower_limit_roi);
    %         message = ['The lower limit criterion for the DTI ROI: ' lower_limit_str ', i.e. ' lower_limit_roi_str ' * std(manually selected area)'];
    %         save_text = [save_text message];
    %         results = [results lower_limit 'lower limit criterion for the DTI ROI' lower_limit_roi '*std(manually selected area)'];
    %         disp(message);



    %% ADC & FA analysis with DTI ROI
    %         ADC_mean_core_str = num2str(ADC_mean_core);
    %         ADC_mean_rim_str = num2str(ADC_mean_rim);
    %         message = ['ADC core and outer part mean value: ' ,ADC_mean_core_str];
    %         save_text = [save_text message];
    %         disp(message);
    %         results = [results ADC_mean_core 'ADC core and outer part mean value' ADC_mean_rim_str 'ADC rim mean value'];
    %         message = ['ADC rim mean value: ', ADC_mean_rim_str];
    %         save_text = [save_text message];
    %         disp(message);
    %         disp('(So far, the core and rim of the ADC image are determined by simply seperating areas of pixel')
    %         disp('with values larger than the mean value of the whole ROI and pixels with values smaller than')
    %         disp('the mean value of the whole ROI.');
    %         FA_mean_core_str = num2str(FA_mean_core);
    %         FA_mean_rim_str = num2str(FA_mean_rim);
    %         message = ['FA core and outer part mean value: ' ,FA_mean_core_str];
    %         save_text = [save_text message];
    %         results = [results FA_mean_core 'FA core and outer part mean value' FA_mean_rim 'FA rim mean value'];
    %         disp(message);
    %         message = ['FA rim mean value: ', FA_mean_rim_str];
    %         save_text = [save_text message];
    %         disp(message);
    %% END ADC & FA analysis with DTI ROI

    my_pixels = [1 2 3 4 5 6];
%     figure;
    %         subplot(4,2,1);plot(my_pixels,core_avg, 'o'); set(gca,'XTick',1:1:6); title('core unshifted');
    %         subplot(4,2,2);plot(my_pixels, rim_avg, 'o'); set(gca,'XTick',1:1:6); title('rim unshifted');
    %         %subplot(5,2,3);plot(my_pixels,core_o_avg, 'o'); set(gca,'XTick',1:1:6); title('core o avg');
    %         %subplot(5,2,4);plot(my_pixels, rim_o_avg, 'o'); set(gca,'XTick',1:1:6); title('rim o avg');
    %         subplot(4,2,3);plot(my_pixels, core_b_avg,'o'); set(gca,'XTick',1:1:6); title('core shifted');
    %         subplot(4,2,4);plot(my_pixels, rim_b_avg, 'o'); set(gca,'XTick',1:1:6); title('rim shifted');
    %         subplot(4,2,5);plot(my_pixels, core_b_smaller_avg,'o'); set(gca,'XTick',1:1:6); title('smaller core shifted');
    %         subplot(4,2,6);plot(my_pixels, rim_b_smaller_avg, 'o'); set(gca,'XTick',1:1:6); title('larger rim shifted');
    %         subplot(4,2,7);plot(my_pixels, tissue_avg,'o'); set(gca,'XTick',1:1:6); title('tissue avg');
    %         subplot(4,2,8);plot(my_pixels, lesion_avg,'o'); set(gca,'XTick',1:1:6); title('lesion avg');

    %         subplot(1,2,1);plot(my_pixels, tissue_avg,'o'); set(gca,'XTick',1:1:6); title('tissue avg');
    %         subplot(1,2,2);plot(my_pixels, lesion_avg,'o'); set(gca,'XTick',1:1:6); title('lesion avg');

%     plot(my_pixels, tissue_avg,'x', my_pixels, tissue2_avg,'+',my_pixels, lesion_avg, 'o'); set(gca,'XTick',1:1:6); title('peripheral areas and lesion average, relative change');
%     h = legend('tissue average','tissue 2 average','lesion average','Location', 'Best'); set(h,'Interpreter','none')
%     savename = [save_folder '/time_series.png'];
%     saveas(gcf,savename)
%     savename = [save_folder '/time_series.fig'];
%     saveas(gcf,savename)
% 
%     plot(my_pixels, d_tissue_avg,'x', my_pixels, d_tissue2_avg,'+',my_pixels, d_lesion_avg, 'o'); set(gca,'XTick',1:1:6); title('peripheral areas and lesion average, difference image');
%     h = legend('tissue average','tissue 2 average','lesion average','Location', 'Best'); set(h,'Interpreter','none')
%     savename = [save_folder '/time_series_diff.png'];
%     saveas(gcf,savename)
%     savename = [save_folder '/time_series_diff.fig'];
%     saveas(gcf,savename)
% 
%     plot(my_pixels, wd_tissue_avg,'x', my_pixels, wd_tissue2_avg,'+',my_pixels, wd_lesion_avg, 'o'); set(gca,'XTick',1:1:6); title('peripheral areas and lesion average, weighted difference');
%     h = legend('tissue average','tissue 2 average','lesion average','Location', 'Best'); set(h,'Interpreter','none')
%     savename = [save_folder '/time_series_wdiff.png'];
%     saveas(gcf,savename)
%     savename = [save_folder '/time_series_wdiff.fig'];
%     saveas(gcf,savename)
% 
%     plot(my_pixels, r_tissue_avg,'x', my_pixels, r_tissue2_avg,'+',my_pixels, r_lesion_avg, 'o'); set(gca,'XTick',1:1:6); title('DCE series mean values');
%     h = legend('tissue average','tissue 2 average','lesion average','Location', 'Best'); set(h,'Interpreter','none')
%     savename = [save_folder '/time_series_basic.png'];
%     saveas(gcf,savename)
%     savename = [save_folder '/time_series_basic.fig'];
%     saveas(gcf,savename)

    %
    %         savename = [save_folder2 '/time_series.png'];
    %         saveas(gcf,savename)
    %         savename = [save_folder2 '/time_series.fig'];
    %         saveas(gcf,savename)
    %my_time = clock;
    %name  = ['Research_results', num2str(my_time(2)),'_',num2str(my_time(3)),'_',num2str(my_time(1)),'_',num2str(my_time(4)),'_',num2str(my_time(5))];

    shift_question =1;
%     if shift_question == 1
%         %% ADC image analysis start
%         shift_test = 0;
% %         figure;
%         while shift_test == 0
% 
% 
%             if shift_x <0
%                 direc_x = -2;
%             else
%                 direc_x = 2;
%             end
%             if shift_y <0
%                 direc_y = -1;
%             else
%                 direc_y = 1;
%             end
%             shift_x_abs = abs(shift_x);
%             shift_y_abs = abs(shift_y);
%             ADC_temp = shift_image(shift_image(ADCimage,shift_y_abs,direc_y),shift_x_abs,direc_x);
%             ADC_temp_small = ADC_temp((y1:y1+y2-1),x1:(x1+x2-1));
%             FA_temp =  shift_image(shift_image(FAimage ,shift_y_abs,direc_y),shift_x_abs,direc_x);
%             FA_temp_small = FA_temp((y1:y1+y2-1),x1:(x1+x2-1));
%             unw_DTI_temp = shift_image(shift_image(unw_DTI,shift_y_abs,direc_y),shift_x_abs,direc_x);
%             unw_DTI_temp = unw_DTI_temp((y1:y1+y2-1),x1:(x1+x2-1));
%             first_DTI_temp = shift_image(shift_image(first_DTI,shift_y_abs,direc_y),shift_x_abs,direc_x);
%             first_DTI_temp = first_DTI_temp((y1:y1+y2-1),x1:(x1+x2-1));
%             avg_DTI_temp = shift_image(shift_image(avg_DTI,shift_y_abs,direc_y),shift_x_abs,direc_x);
%             avg_DTI_temp = avg_DTI_temp((y1:y1+y2-1),x1:(x1+x2-1));
%             %                 ADC_temp_DCE_ROI = ADC_temp.*sur_and_mar_inv;
% %             subplot(5,2,1);imshow(ADC_temp_DCE_ROI.*sur_and_mar_inv, 'InitialMagnification', 'fit'); colorbar;  caxis([min(ADC_temp_DCE_ROI(:)) max(ADC_temp_DCE_ROI(:))]); %colormap(jet)
% %             subplot(5,2,2);imshow(ADC_temp_small.*sur_and_mar_inv, 'InitialMagnification', 'fit'); colorbar; caxis([min(ADC_temp(:)) max(ADC_temp(:))]);
% %             subplot(5,2,3);imshow(FA_temp_DCE_ROI.*sur_and_mar_inv, 'InitialMagnification', 'fit'); colorbar; caxis([min(FA_temp_DCE_ROI(:)) max(FA_temp_DCE_ROI(:))]);
% %             subplot(5,2,4);imshow(FA_temp_small.*sur_and_mar_inv, 'InitialMagnification', 'fit'); colorbar; caxis([min(FA_temp(:)) max(FA_temp(:))]);
% %             subplot(5,2,5);imshow(unw_DTI_temp_DCE_ROI.*sur_and_mar_inv, 'InitialMagnification', 'fit'); colorbar; caxis([min(unw_DTI_temp_DCE_ROI(:)) max(unw_DTI_temp_DCE_ROI(:))]);
% %             subplot(5,2,6);imshow(unw_DTI_temp.*sur_and_mar_inv, 'InitialMagnification', 'fit'); colorbar; caxis([min(unw_DTI_temp(:)) max(unw_DTI_temp(:))]);
% %             subplot(5,2,7);imshow(first_DTI_temp_DCE_ROI.*sur_and_mar_inv, 'InitialMagnification', 'fit'); colorbar; caxis([min(first_DTI_temp_DCE_ROI(:)) max(first_DTI_temp_DCE_ROI(:))]);
% %             subplot(5,2,8);imshow(first_DTI_temp.*sur_and_mar_inv, 'InitialMagnification', 'fit'); colorbar; caxis([min(first_DTI_temp(:)) max(first_DTI_temp(:))]);
% %             subplot(5,2,9);imshow(avg_DTI_temp_DCE_ROI.*sur_and_mar_inv, 'InitialMagnification', 'fit'); colorbar; caxis([min(avg_DTI_temp_DCE_ROI(:)) max(avg_DTI_temp_DCE_ROI(:))]);
% %             subplot(5,2,10);imshow(avg_DTI_temp.*sur_and_mar_inv, 'InitialMagnification', 'fit'); colorbar; caxis([min(avg_DTI_temp(:)) max(avg_DTI_temp(:))]);
% 
%                 shift_test = 1; % the DTI/DCE should be the same for all slices
% 
%         end
% 
%         shift_word1 = num2str(shift_x_abs);
%         if shift_x < 0
%             shift_word1 = ['minus_' shift_word1];
%         end
% 
%         shift_word2 = num2str(shift_y_abs);
%         if shift_y < 0
%             shift_word2 = ['minus_' shift_word2];
%         end
% 
%         savename = [save_folder '/DTI_shift_x_' shift_word1 '_y_' shift_word2 '.fig']; % save as fig...
%         saveas(gcf,savename)
%         savename = [save_folder '/DTI_shift_x_' shift_word1 '_y_' shift_word2 '.png']; % ... and png.
%         saveas(gcf,savename)
% 
%         %             savename = [save_folder2 '/DTI_shift_x_' shift_word1 '_y_' shift_word2 '.fig']; % save as fig...
%         %             saveas(gcf,savename)
%         %             savename = [save_folder2 '/DTI_shift_x_' shift_word1 '_y_' shift_word2 '.png']; % ... and png.
%         %             saveas(gcf,savename)
%         %
% %         figure; imshow(ADC_temp_small,'InitialMagnification', 'fit'); colorbar; caxis([0 3]);
% % 
% %         savename = [save_folder '/ADC_with_DCE_ROI_shifted.fig']; % save as fig...
% %         saveas(gcf,savename)
% %         savename = [save_folder '/ADC_with_DCE_ROI_shifted.png']; % ... and png.
% %         saveas(gcf,savename)
%         %
%         %             savename = [save_folder2 '/ADC_with_DCE_ROI_shifted.fig']; % save as fig...
%         %             saveas(gcf,savename)
%         %             savename = [save_folder2 '/ADC_with_DCE_ROI_shifted.png']; % ... and png.
%         %             saveas(gcf,savename)
% 
%         %% ADC analysis with DCE ROI
% 
% 
%         ADC_tissue_shifted = ADC_temp.*DCE_tissue;
%         ADC_tissue_vec_shifted = [ADC_tissue_vec_shifted ADC_tissue_shifted(:)'];
%         if slice == nslices
%             ADC_tissue_mean_shifted = mean(ADC_tissue_vec_shifted(ADC_tissue_vec_shifted~=0));
%             ADC_tissue_std_shifted = std(ADC_tissue_vec_shifted(ADC_tissue_vec_shifted~=0));
%         end
% 
%         ADC_tissue2_shifted = ADC_temp.*DCE_tissue2;
%         ADC_tissue2_vec_shifted = [ADC_tissue2_vec_shifted ADC_tissue2_shifted(:)'];
%         if slice == nslices
%             ADC_tissue2_mean_shifted = mean(ADC_tissue2_vec_shifted(ADC_tissue2_vec_shifted~=0));
%             ADC_tissue2_std_shifted = std(ADC_tissue2_vec_shifted(ADC_tissue2_vec_shifted~=0));
%         end
% 
%         ADC_lesion_shifted = ADC_temp.*DCE_inner;
%         ADC_lesion_vec_shifted = [ADC_lesion_vec_shifted ADC_lesion_shifted(:)'];
%         if slice == nslices
%             ADC_lesion_mean_shifted = mean(ADC_lesion_vec_shifted(ADC_lesion_vec_shifted~=0));
%             ADC_lesion_std_shifted = std(ADC_lesion_vec_shifted(ADC_lesion_vec_shifted~=0));
%         end
% 
%         %% END ADC analysis with DCE ROI
% 
%         %%% FA image analysis start
%         %         [filename,pathname] = uigetfile({'*.m;*.mat','MATLAB Files (*.m,*.mat)'},'Select the slice corresponding FA image','/home/hahntobi/matlab/patient_images/');
%         %         file = [pathname filename];
% 
% 
% %         figure; imshow(FA_temp_small.*sur_and_mar_inv, 'InitialMagnification', 'fit'); colorbar; caxis([0 1.2]); %colormap(jet)
% % 
% %         savename = [save_folder '/FA_with_DCE_ROI_shifted.fig']; % save as fig...
% %         saveas(gcf,savename)
% %         savename = [save_folder '/FA_with_DCE_ROI_shifted.png']; % ... and png.
% %         saveas(gcf,savename)
% 
%         %             savename = [save_folder2 '/FA_with_DCE_ROI_shifted.fig']; % save as fig...
%         %             saveas(gcf,savename)
%         %             savename = [save_folder2 '/FA_with_DCE_ROI_shifted.png']; % ... and png.
%         %             saveas(gcf,savename)
% 
%         FA_tissue_shifted = FA_temp.*DCE_tissue;
%         FA_tissue_vec_shifted = [FA_tissue_vec_shifted FA_tissue_shifted(:)'];
%         if slice == nslices
%             FA_tissue_mean_shifted = mean(FA_tissue_vec_shifted(FA_tissue_vec_shifted~=0));
%             FA_tissue_std_shifted = std(FA_tissue_vec_shifted(FA_tissue_vec_shifted~=0));
%         end
% 
%         FA_tissue2_shifted = FA_temp.*DCE_tissue2;
%         FA_tissue2_vec_shifted = [FA_tissue2_vec_shifted FA_tissue2_shifted(:)'];
%         if slice == nslices
%             FA_tissue2_mean_shifted = mean(FA_tissue2_vec_shifted(FA_tissue2_vec_shifted~=0));
%             FA_tissue2_std_shifted = std(FA_tissue2_vec_shifted(FA_tissue2_vec_shifted~=0));
%         end
% 
%         FA_lesion_shifted = FA_temp.*DCE_inner;
%         FA_lesion_vec_shifted = [FA_lesion_vec_shifted FA_lesion_shifted(:)'];
%         if slice == nslices
%             FA_lesion_mean_shifted = mean(FA_lesion_vec_shifted(FA_lesion_vec_shifted~=0));
%             FA_lesion_std_shifted = std(FA_lesion_vec_shifted(FA_lesion_vec_shifted~=0));
%         end
%     end
% 
end
% %% ADC & FA analysis with DCE ROI
% %results = [results FA_ROI_DCE_mean 'FA_ROI_DCE_mean' FA_ROI_DCE_std 'FA_ROI_DCE_std' ADC_ROI_DCE_mean 'ADC_ROI_DCE_mean' ADC_ROI_DCE_std 'ADC_ROI_DCE_std'];
% %FA_ROI_DCE_mean_str = num2str(FA_ROI_DCE_mean);
% %FA_ROI_DCE_std_str = num2str(FA_ROI_DCE_std);
% %ADC_ROI_DCE_mean_str = num2str(ADC_ROI_DCE_mean);
% %ADC_ROI_DCE_std_str = num2str(ADC_ROI_DCE_std);
% %message = ['FA image with DCE ROI: ' FA_ROI_DCE_mean_str ' +- ' FA_ROI_DCE_std_str '. ADC image with DCE ROI: ' ADC_ROI_DCE_mean_str ' +- ' ADC_ROI_DCE_std_str];
% %disp(message);
% results = [results FA_tissue_mean 'FA tissue mean' FA_tissue_std 'FA tissue std' FA_tissue2_mean 'FA tissue 2 mean' FA_tissue2_std 'FA tissue 2 std' FA_lesion_mean 'FA lesion mean' FA_lesion_std 'FA lesion std' ];
% results = [results ADC_tissue_mean 'ADC tissue mean' ADC_tissue_std 'ADC tissue std' ADC_tissue2_mean 'ADC tissue 2 mean' ADC_tissue2_std 'ADC tissue 2 std' ADC_lesion_mean 'ADC lesion mean' ADC_lesion_std 'ADC lesion std'];
% %% END ADC & FA analysis with DCE ROI

%% DTI analysis with DCE ROI
%results = [results DTI_unw_ROI_mean 'DTI_unweighted_DCE_mean' DTI_first_ROI_mean 'DTI_first_DCE_mean'];
%% END DTI analysis with DCE ROI

%% WASH-IN & WASH-OUT analysis with DCE ROI
results = [results wash_in_ROI_DCE_mean 'wash-in lesion mean' wash_in_ROI_DCE_std 'wash-in lesion std' wash_in_tissue_mean 'wash-in tissue mean' wash_in_tissue_std 'wash-in tissue std' wash_in_tissue2_mean 'wash-in tissue 2 mean' wash_in_tissue2_std 'wash-in tissue 2 std'];
results = [results wash_out_ROI_DCE_mean 'wash-out lesion mean' wash_out_ROI_DCE_std 'wash-out lesion std' wash_out_tissue_mean 'wash-out tissue mean' wash_out_tissue_std 'wash-out tissue std' wash_out_tissue2_mean 'wash-out tissue 2 mean' wash_out_tissue2_std 'wash-out tissue 2 std'];
results = [results wash_out_short_ROI_DCE_mean 'wash-out_short lesion mean' wash_out_short_ROI_DCE_std 'wash-out_short lesion std' wash_out_short_tissue_mean 'wash-out_short tissue mean' wash_out_short_tissue_std 'wash-out_short tissue std' wash_out_short_tissue2_mean 'wash-out_short tissue 2 mean' wash_out_short_tissue2_std 'wash-out_short tissue 2 std'];
results = [results DCEimage_lesion_mean(1) 'DCE image no.1 lesion mean' DCEimage_lesion_mean(2) 'DCE image no.2 lesion mean' DCEimage_lesion_mean(3) 'DCE image no.3 lesion mean' DCEimage_lesion_mean(4) 'DCE image no.4 lesion mean' DCEimage_lesion_mean(5) 'DCE image no.5 lesion mean' DCEimage_lesion_mean(6) 'DCE image no.6 lesion mean'];
results = [results DCEimage_tissue_mean(1) 'DCE image no.1 tissue mean' DCEimage_tissue_mean(2) 'DCE image no.2 tissue mean' DCEimage_tissue_mean(3) 'DCE image no.3 tissue mean' DCEimage_tissue_mean(4) 'DCE image no.4 tissue mean' DCEimage_tissue_mean(5) 'DCE image no.5 tissue mean' DCEimage_tissue_mean(6) 'DCE image no.6 tissue mean'];
results = [results DCEimage_tissue2_mean(1) 'DCE image no.1 tissue2 mean' DCEimage_tissue2_mean(2) 'DCE image no.2 tissue2 mean' DCEimage_tissue2_mean(3) 'DCE image no.3 tissue2 mean' DCEimage_tissue2_mean(4) 'DCE image no.4 tissue2 mean' DCEimage_tissue2_mean(5) 'DCE image no.5 tissue2 mean' DCEimage_tissue2_mean(6) 'DCE image no.6 tissue2 mean'];
results = [results DCEimage_lesion_std(1) 'DCE image no.1 lesion std' DCEimage_lesion_std(2) 'DCE image no.2 lesion std' DCEimage_lesion_std(3) 'DCE image no.3 lesion std' DCEimage_lesion_std(4) 'DCE image no.4 lesion std' DCEimage_lesion_std(5) 'DCE image no.5 lesion std' DCEimage_lesion_std(6) 'DCE image no.6 lesion std'];
results = [results DCEimage_tissue_std(1) 'DCE image no.1 tissue std' DCEimage_tissue_std(2) 'DCE image no.2 tissue std' DCEimage_tissue_std(3) 'DCE image no.3 tissue std' DCEimage_tissue_std(4) 'DCE image no.4 tissue std' DCEimage_tissue_std(5) 'DCE image no.5 tissue std' DCEimage_tissue_std(6) 'DCE image no.6 tissue std'];
results = [results DCEimage_tissue2_std(1) 'DCE image no.1 tissue2 std' DCEimage_tissue2_std(2) 'DCE image no.2 tissue2 std' DCEimage_tissue2_std(3) 'DCE image no.3 tissue2 std' DCEimage_tissue2_std(4) 'DCE image no.4 tissue2 std' DCEimage_tissue2_std(5) 'DCE image no.5 tissue2 std' DCEimage_tissue2_std(6) 'DCE image no.6 tissue2 std'];
results = [results diff_lesion_mean 'Difference image lesion mean' diff_tissue_mean 'Difference image tissue mean' diff_tissue2_mean 'Difference image tissue2 mean'];
results = [results diff_lesion_std 'Difference image lesion std' diff_tissue_std 'Difference image tissue std' diff_tissue2_std 'Difference image tissue2 std'];
results = [results wdiff_lesion_mean 'Weighted difference image lesion mean' wdiff_tissue_mean 'Weighted difference image tissue mean' wdiff_tissue2_mean 'Weighted difference image tissue2 mean'];
results = [results wdiff_lesion_std 'Weighted difference image lesion std' wdiff_tissue_std 'Weighted difference image tissue std' wdiff_tissue2_std 'Weighted difference image tissue2 std'];
%% END WASH-IN & WASH-OUT analysis with DCE ROI


%% append the optional shifted ADC/FA values
% if shift_question == 1
%     results = [results FA_tissue_mean_shifted 'shifted FA tissue mean' FA_tissue_std_shifted 'shifted FA tissue std' FA_tissue2_mean_shifted 'shifted FA tissue 2 mean' FA_tissue2_std_shifted 'shifted FA tissue 2 std' FA_lesion_mean_shifted 'shifted FA lesion mean' FA_lesion_std_shifted 'shifted FA lesion std' ];
%     results = [results ADC_tissue_mean_shifted 'shifted ADC tissue mean' ADC_tissue_std_shifted 'shifted ADC tissue std' ADC_tissue2_mean_shifted 'shifted ADC tissue 2 mean' ADC_tissue2_std_shifted 'shifted ADC tissue 2 std' ADC_lesion_mean_shifted 'shifted ADC lesion mean' ADC_lesion_std_shifted 'shifted ADC lesion std' ];
% end
%% end appended shifted ADC/FA values
results2 = cell(1);
for i =1:((size(results,2)-1)/2); results2(i,1)=results(1,2*i);end
for i =1:((size(results,2)-1)/2); results2(i,2)=results(1,1+2*i);end
results = results2;

results_copy = results;

savename = [save_folder '/analysis_results.mat'];
save(savename, 'results');
%         savename = [save_folder2 '/analysis_results.mat'];
%         save(savename, 'results');
copy_folder = ['/home/hahntobi/matlab/patient_images/' SID '/'];
copyfile(copy_folder,save_folder); % copy all images
%         copyfile(copy_folder,save_folder2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DETAILED ANALYSIS


% open the specific folder
dir_name = save_folder_basic;
% temp3 = [dir_name '/research_results/slice3'];
% temp2 = [dir_name '/research_results/slice2'];
% temp1 = [dir_name '/research_results/slice1'];
% flag1 = 0;
% flag2 = 0;
% flag3 = 0;
% 
% if exist(temp3)
%     dir_data = temp3;
%     flag3 = 1;
%     temp = [dir_data '/final_BWs.mat'];
%     load(temp);
%     inner_BW_large3 = inner_BW_large;
% elseif exist(temp2)
%     dir_data = temp2;
% elseif exist(temp1)
%     dir_data = temp1;
% end
% 
% if exist(temp2)
%     flag2 = 1;
%     temp = [temp2 '/final_BWs.mat'];
%     load(temp);
%     inner_BW_large2 = inner_BW_large;
% end
% if exist(temp1)
%     flag1 = 1;
%     temp = [temp1 '/final_BWs.mat'];
%     load(temp);
%     inner_BW_large1 = inner_BW_large;
% end

figure
rec = gcf;

% display('new patient');
mal = zeros(3,1);
ben = zeros(3,1);
plat= zeros(3,1);
% 
% shift_x = input('Please indicate the shift in the positive x-direction (Matlab convention): ');
% shift_y = input('Please indicate the shift in the positive y-direction (Matlab convention): ');

if shift_x <0
    direc_x = -2;
else
    direc_x = 2;
end
if shift_y <0
    direc_y = -1;
else
    direc_y = 1;
end
shift_x_abs = abs(shift_x);
shift_y_abs = abs(shift_y);

%% THE CUT-OFF VALUES

col = -15; % cut-off value on the "left side" (standard: -30 and 30, resp.)
cor = +30; % cut-off value on the "right side"

%%

mal_mean = []; % '3D analysis'
ben_mean = [];
plat_mean = [];

BW_wo_vol = 0;
BW_plat_vol = 0;
BW_pers_vol = 0;

% if flag1 ==1
%     temp = fullfile(dir_data, '/DCE_temp/wash_out1.mat');
%     wash_out1 = load(temp, 'image');
%     wash_out1 = wash_out1.image;
%     wash_out1 = wash_out1.*inner_BW_large1;
%     BW_wo1 = roicolor(wash_out1, -90,col).*inner_BW_large1;
%     BW_plat1 = roicolor(wash_out1, col, cor).*inner_BW_large1;
%     BW_pers1 = roicolor(wash_out1, cor, 90).*inner_BW_large1;
%     BW_wo1_area = sum(sum(BW_wo1==1));
%     BW_plat1_area = sum(sum(BW_plat1==1));
%     BW_pers1_area = sum(sum(BW_pers1==1));
%     BW_wo_vol = BW_wo_vol + BW_wo1_area;
%     BW_plat_vol = BW_plat_vol + BW_plat1_area;
%     BW_pers_vol = BW_pers_vol + BW_pers1_area;
% 
% 
% end
% 
% if flag2 ==1
%     temp = fullfile(dir_data, '/DCE_temp/wash_out2.mat');
%     wash_out2 = load(temp, 'image');
%     wash_out2 = wash_out2.image;
%     wash_out2 = wash_out2.*inner_BW_large2;
%     BW_wo2 = roicolor(wash_out2, -90,col).*inner_BW_large2;
%     BW_pers2 = roicolor(wash_out2, cor, 90).*inner_BW_large2;
%     BW_plat2 = roicolor(wash_out2, col, cor).*inner_BW_large2;
%     BW_wo2_area = sum(sum(BW_wo2==1));
%     BW_plat2_area = sum(sum(BW_plat2==1));
%     BW_pers2_area = sum(sum(BW_pers2==1));  
%         BW_wo_vol = BW_wo_vol + BW_wo2_area;
%     BW_plat_vol = BW_plat_vol + BW_plat2_area;
%     BW_pers_vol = BW_pers_vol + BW_pers2_area;
% 
% end
% 
% if flag3 ==1
%     temp = fullfile(dir_data, '/DCE_temp/wash_out3.mat');
%     wash_out3 = load(temp, 'image');
%     wash_out3 = wash_out3.image;
%     wash_out3 = wash_out3.*inner_BW_large3;
%     BW_wo3 = roicolor(wash_out3, -90,col).*inner_BW_large3;
%     BW_pers3 = roicolor(wash_out3, cor, 90).*inner_BW_large3;
%     BW_plat3 = roicolor(wash_out3, col, cor).*inner_BW_large3;
%     BW_wo3_area = sum(sum(BW_wo3==1));
%     BW_plat3_area = sum(sum(BW_plat3==1));
%     BW_pers3_area = sum(sum(BW_pers3==1));    
%         BW_wo_vol = BW_wo_vol + BW_wo3_area;
%     BW_plat_vol = BW_plat_vol + BW_plat3_area;
%     BW_pers_vol = BW_pers_vol + BW_pers3_area;
% 
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FA_ROI_vec = [];
ADC_ROI_vec = [];
ADC_pers_vec = [];
ADC_plat_vec = [];
ADC_wo_vec = [];
FA_pers_vec = [];
FA_plat_vec = [];
FA_wo_vec = [];
ADC_pers_vec_shifted = [];
ADC_plat_vec_shifted = [];
ADC_wo_vec_shifted = [];
FA_pers_vec_shifted = [];
FA_plat_vec_shifted = [];
FA_wo_vec_shifted = [];
DTI_unw_ROI_vec = [];
DTI_first_ROI_vec = [];
wash_in_wo_vec = [];
wash_in_pers_vec = [];
wash_in_plat_vec = [];
wash_out_wo_vec = [];
wash_out_pers_vec = [];
wash_out_plat_vec = [];
diff_wo_vec = [];
diff_pers_vec = [];
diff_plat_vec = [];
wdiff_wo_vec = [];
wdiff_pers_vec = [];
wdiff_plat_vec = [];
wash_out_short_ROI_vec = [];
wash_out_short_pers_vec = [];
wash_out_short_plat_vec = [];
DCEimage_pers_vec = zeros(6,3*512*512);
DCEimage_plat_vec = zeros(6,3*512*512);
DCEimage_wo_vec = zeros(6,3*512*512);
results = cell(1); % result - cell

BW_wo_volume = 0;
BW_plat_volume = 0;
BW_pers_volume = 0;

%%% SAVE FOLDER
save_folder_basic = dir_name;
% save_folder_basic = uigetdir('/export/res/breast','Please select the
% patient subject folder or the main save folder for the analysis results');
%save_folder2_basic = uigetdir('/home/hahntobi/matlab/results/', 'Please select a copy save folder outside of vali');
%%% END SAVE FOLDER


for slice = 1:nslices

    if slice == 1
        user_entry_slice = user_entry_slice_copy;
    else
        user_entry_slice = user_entry_slice+1; % The case that only slice 1 and slice 3 are examined shouldn't happen!
    end
    user_entry_slice_str = num2str(user_entry_slice);
    slicen =num2str(slice);
    save_folder = [save_folder_basic '/research_results_detailed/slice' user_entry_slice_str];
    exist_test = exist(save_folder,'dir');
    if exist_test == 0
        mkdir(save_folder);
    end

    s1 = ['/research_results_detailed/slice' user_entry_slice_str];
    temp_dir = [dir_name s1];
    
    last_slice = first_slice+nslices-1;
    last_slice_str = num2str(last_slice);
    s1 = ['/research_results_detailed/slice' last_slice_str];
    dir_data = [dir_name s1];
    
%     s1 =['temp_dir = [dir_name "/research_results/slice" user_entry_slice_str];'];
%     eval(s1);
%     s2 = ['temp = [temp1 "/final_BWs.mat"]'];
%     eval(s2);
    temp = [temp_dir '/final_BWs.mat'];
    load(temp);
    
    temp = ['/DCE_temp/wash_out' user_entry_slice_str '.mat'];
    temp = fullfile(dir_data, temp);

%     s2 = ['temp = fullfile(dir_data, "/DCE_temp/wash_out" user_entry_slice_str ".mat");'];
%     eval(s)

%     temp = fullfile(dir_data, '/DCE_temp/wash_out3.mat');
    wash_out = load(temp, 'image');
    wash_out = wash_out.image;
    wash_out = wash_out.*inner_BW_large;
    BW_wo = roicolor(wash_out, -90,col).*inner_BW_large;
    BW_pers = roicolor(wash_out, cor, 90).*inner_BW_large;
    BW_plat = roicolor(wash_out, col, cor).*inner_BW_large;
    BW_wo_area = sum(sum(BW_wo==1));
    BW_plat_area = sum(sum(BW_plat==1));
    BW_pers_area = sum(sum(BW_pers==1));



    %     if user_entry_slice == 1
    %         BW_wo = BW_wo1;
%         BW_pers = BW_pers1;
%         BW_plat = BW_plat1;
%         
%     elseif user_entry_slice == 2 
%         BW_wo = BW_wo2;
%         BW_pers = BW_pers2;
%         BW_plat = BW_plat2;
% 
%     elseif user_entry_slice == 3
%         BW_wo = BW_wo3;
%         BW_pers = BW_pers3;
%         BW_plat = BW_plat3;
% 
%     end
    
    DCE_image_dir = ['/home/hahntobi/matlab/patient_images/' SID '/DCE_temp'];
    for i=1:6 % six DCE images
        DCE_number = 6*(user_entry_slice - 1) + i;
        DCE_number = num2str(DCE_number);
        if (6*(user_entry_slice-1)+i)<10
            DCE_number = ['0' DCE_number];
        end
        file = [DCE_image_dir '/' DCE_number];
        load(file,'image');
        DCE(:,:,i) = image;
    end

    DCE_ref = double(DCE(:,:,1));
    DCE_ref(DCE_ref==0)=1; % setting zeros to one, we will forget these pixels later by using the 'no_zero_mask', just set so that the division by DCE_ref is possible
    for i=1:6
        DCE_temp = double(DCE(:,:,i)); % double-conversion to make sure that the calculation of the std is correct
        

        %% RAW ('BASIC') VALUES CALULATION FOR TIME_SERIES (just for
        %% plotting)
        mask = BW_pers;
        temp = DCE_temp.*mask;
        temp_vec = temp(mask~=0);
        temp_vec(temp_vec==0)=[]; % regard the case that DCE_temp == 0
        r_pers_avg(i)=mean(temp_vec);
        r_pers_std(i)=std(temp_vec);
        mask = BW_wo;
        temp = DCE_temp.*mask;
        temp_vec = temp(mask~=0);
        temp_vec(temp_vec==0)=[];
        r_wo_avg(i)=mean(temp_vec);
        r_wo_std(i)=std(temp_vec);
        mask = BW_plat;
        temp = DCE_temp.*mask;
        temp_vec = temp(mask~=0);
        temp_vec(temp_vec==0)=[];
        r_plat_avg(i)=mean(temp_vec);
        r_plat_std(i)=std(temp_vec);    

        %% END RAW VALUES CALCULATION FOR TIME_SERIES


        DCEimage_pers = DCE_temp.*BW_pers;
        DCEimage_pers_vec(i,(user_entry_slice-1)*512*512+1:user_entry_slice*512*512) = DCEimage_pers(:)';
        if slice == nslices
            temp = DCEimage_pers_vec(i,:);
            DCEimage_pers_mean(i) = mean(temp(temp~=0));
            DCEimage_pers_std(i) = std(temp(temp~=0));
        end

        DCEimage_plat = DCE_temp.*BW_plat;
        DCEimage_plat_vec(i,(user_entry_slice-1)*512*512+1:user_entry_slice*512*512) = DCEimage_plat(:)';
        if slice == nslices
            temp = DCEimage_plat_vec(i,:);
            DCEimage_plat_mean(i) = mean(temp(temp~=0));
            DCEimage_plat_std(i) = std(temp(temp~=0));
        end

        DCEimage_wo = DCE_temp.*BW_wo;
        DCEimage_wo_vec(i,(user_entry_slice-1)*512*512+1:user_entry_slice*512*512) = DCEimage_wo(:)';
        if slice == nslices
            temp = DCEimage_wo_vec(i,:);
            DCEimage_wo_mean(i) = mean(temp(temp~=0));
            DCEimage_wo_std(i) = std(temp(temp~=0));
        end
        %         figure; subplot(3,2,1); imshow(DCE(:,:,i).*DTI_inner_o_shift.*BW_wo,[]);
        %         subplot(3,2,2); imshow(DCE(:,:,i).*(BW_wo-DTI_inner_o_shift.*BW_wo),[]);
        %         subplot(3,2,3); imshow(DCE(:,:,i).*DTI_inner_b_shift.*BW_wo,[]);
        %         subplot(3,2,4); imshow(DCE(:,:,i).*(BW_wo-DTI_inner_b_shift.*BW_wo),[]);
        %         subplot(3,2,5); imshow(DCE(:,:,i).*BW_pers,[]);
    end

    %         results = [results core_b_avg 'core_b_avg' core_std 'core_std' rim_b_avg 'rim_b_avg' rim_std 'rim_std' core_b_smaller_avg 'core_b_smaller_avg' core_b_smaller_std 'core_b_smaller_std' rim_b_smaller_avg 'rim_b_smaller_avg' rim_b_smaller_std 'rim_b_smaller_std' pers_avg 'pers_avg' pers_std 'pers_std' wo_avg 'wo_avg' wo_std 'wo_std'];

    %%% ADC image analysis start
%     file = ['/home/hahntobi/matlab/patient_images/' SID '/ADC_temp/0' user_entry_slice_str '.mat'];
%     load(file,'image');
%     ADCimage = image;
%     ADCimage = size_adapt(ADCimage);
%     ADCimage = shift_image(shift_image(ADCimage,shift_y_abs,direc_y),shift_x_abs,direc_x);
% 
%     ADC_pers = ADCimage.*BW_pers;
%     ADC_pers_vec = [ADC_pers_vec ADC_pers(:)'];
%     if slice == nslices
%         ADC_pers_mean = mean(ADC_pers_vec(ADC_pers_vec~=0));
%         ADC_pers_std = std(ADC_pers_vec(ADC_pers_vec~=0));
%     end
% 
%     ADC_plat = ADCimage.*BW_plat;
%     ADC_plat_vec = [ADC_plat_vec ADC_plat(:)'];
%     if slice == nslices
%         ADC_plat_mean = mean(ADC_plat_vec(ADC_plat_vec~=0));
%         ADC_plat_std = std(ADC_plat_vec(ADC_plat_vec~=0));
%     end
% 
%     ADC_wo = ADCimage.*BW_wo;
%     ADC_wo_vec = [ADC_wo_vec ADC_wo(:)'];
%     if slice == nslices
%         ADC_wo_mean = mean(ADC_wo_vec(ADC_wo_vec~=0));
%         ADC_wo_std = std(ADC_wo_vec(ADC_wo_vec~=0));
%     end
% 
%     %% END ADC analysis with DCE ROI
% 
% 
%     %%% FA image analysis start
%     file = ['/home/hahntobi/matlab/patient_images/' SID '/FA_temp/0' user_entry_slice_str '.mat'];
%     load(file,'image');
%     FAimage = image;
%     FAimage = size_adapt(FAimage);
%     FAimage = shift_image(shift_image(FAimage,shift_y_abs,direc_y),shift_x_abs,direc_x);
%     
%     FA_pers = FAimage.*BW_pers;
%     FA_pers_vec = [FA_pers_vec FA_pers(:)'];
%     if slice == nslices
%         FA_pers_mean = mean(FA_pers_vec(FA_pers_vec~=0));
%         FA_pers_std = std(FA_pers_vec(FA_pers_vec~=0));
%     end
% 
%     FA_plat = FAimage.*BW_plat;
%     FA_plat_vec = [FA_plat_vec FA_plat(:)'];
%     if slice == nslices
%         FA_plat_mean = mean(FA_plat_vec(FA_plat_vec~=0));
%         FA_plat_std = std(FA_plat_vec(FA_plat_vec~=0));
%     end
% 
%     FA_wo = FAimage.*BW_wo;
%     FA_wo_vec = [FA_wo_vec FA_wo(:)'];
%     if slice == nslices
%         size(FA_wo_vec)
%         FA_wo_mean = mean(FA_wo_vec(FA_wo_vec~=0));
%         FA_wo_std = std(FA_wo_vec(FA_wo_vec~=0));
%     end
%     %% END FA image analysis



    %% WASH-IN Analysis
    file = ['/home/hahntobi/matlab/patient_images/' SID '/DCE_temp/wash_in_slice' user_entry_slice_str '.mat'];
    load(file,'image');
    wash_in = image;
    wash_in_wo = wash_in.*BW_wo;
    wash_in_wo_single_vec = wash_in_wo(:);

    wash_in_wo_vec = [wash_in_wo_vec wash_in_wo_single_vec'];
    if slice == nslices
        wash_in_wo_mean = mean(wash_in_wo_vec(wash_in_wo_vec~=0));
        wash_in_wo_std = std(wash_in_wo_vec(wash_in_wo_vec~=0));
    end

    wash_in_pers = wash_in.*BW_pers;
    wash_in_pers_vec = [wash_in_pers_vec wash_in_pers(:)'];
    if slice == nslices
        wash_in_pers_mean = mean(wash_in_pers_vec(wash_in_pers_vec~=0));
        wash_in_pers_std = std(wash_in_pers_vec(wash_in_pers_vec~=0));
    end

    wash_in_plat = wash_in.*BW_plat;
    wash_in_plat_vec = [wash_in_plat_vec wash_in_plat(:)'];
    if slice == nslices
        wash_in_plat_mean = mean(wash_in_plat_vec(wash_in_plat_vec~=0));
        wash_in_plat_std = std(wash_in_plat_vec(wash_in_plat_vec~=0));
    end

    %% END WASH-IN Analysis

    %% WASH-OUT Analysis
    file = ['/home/hahntobi/matlab/patient_images/' SID '/DCE_temp/wash_out' user_entry_slice_str '.mat'];
    load(file,'image');
    wash_out = image;
    wash_out_wo = wash_out.*BW_wo;
    wash_out_wo_single_vec = wash_out_wo(:);
    
    wash_out_wo_vec = [wash_out_wo_vec wash_out_wo_single_vec'];
    if slice == nslices
        wash_out_wo_mean = mean(wash_out_wo_vec(wash_out_wo_vec~=0));
        wash_out_wo_std = std(wash_out_wo_vec(wash_out_wo_vec~=0));
    end

    wash_out_pers = wash_out.*BW_pers;
    wash_out_pers_vec = [wash_out_pers_vec wash_out_pers(:)'];
    if slice == nslices
        wash_out_pers_mean = mean(wash_out_pers_vec(wash_out_pers_vec~=0));
        wash_out_pers_std = std(wash_out_pers_vec(wash_out_pers_vec~=0));
    end

    wash_out_plat = wash_out.*BW_plat;
    wash_out_plat_vec = [wash_out_plat_vec wash_out_plat(:)'];
    if slice == nslices
        wash_out_plat_mean = mean(wash_out_plat_vec(wash_out_plat_vec~=0));
        wash_out_plat_std = std(wash_out_plat_vec(wash_out_plat_vec~=0));
    end

    %% END WASH-OUT Analysis

    %% 'SHORT' WASH-OUT Analysis
    file = ['/home/hahntobi/matlab/patient_images/' SID '/DCE_temp/wash_out_short' user_entry_slice_str '.mat'];
    load(file,'image');
    wash_out_short = image;
    wash_out_short_ROI = wash_out_short.*BW_wo;
    wash_out_short_ROI_single_vec = wash_out_short_ROI(:);

    wash_out_short_ROI_vec = [wash_out_short_ROI_vec wash_out_short_ROI_single_vec];
    if slice == nslices
        wash_out_short_ROI_DCE_mean = mean(wash_out_short_ROI_vec(wash_out_short_ROI_vec~=0));
        wash_out_short_ROI_DCE_std = std(wash_out_short_ROI_vec(wash_out_short_ROI_vec~=0));
    end

    wash_out_short_pers = wash_out_short.*BW_pers;
    wash_out_short_pers_vec = [wash_out_short_pers_vec wash_out_short_pers(:)];
    if slice == nslices
        wash_out_short_pers_mean = mean(wash_out_short_pers_vec(wash_out_short_pers_vec~=0));
        wash_out_short_pers_std = std(wash_out_short_pers_vec(wash_out_short_pers_vec~=0));
    end

    wash_out_short_plat = wash_out_short.*BW_plat;
    wash_out_short_plat_vec = [wash_out_short_plat_vec wash_out_short_plat(:)];
    if slice == nslices
        wash_out_short_plat_mean = mean(wash_out_short_plat_vec(wash_out_short_plat_vec~=0));
        wash_out_short_plat_std = std(wash_out_short_plat_vec(wash_out_short_plat_vec~=0));
    end

    %% END 'SHORT' WASH-OUT Analysis


    disp('____________________________________________________________________');
    disp('                      ANALYSIS RESULTS                              '); 

    % kinetic curves plots
    my_pixels = [1 2 3 4 5 6];
    figure;
    plot(my_pixels, r_pers_avg,'x', my_pixels, r_plat_avg,'+',my_pixels, r_wo_avg, 'o'); set(gca,'XTick',1:1:6); title('DCE series mean values for wash out, plateau and persistent areas');
    h = legend('pers average','plateau average','wo average','Location', 'Best'); set(h,'Interpreter','none')
    %%% (if-clause is not really necessary, it's automatically sorted into
    %%% the right slice-directory...)
%     if user_entry_slice == 1
%         savename = [save_folder '/time_series_slice1_detailed_analysis.png'];
%         saveas(gcf,savename)
%         savename = [save_folder '/time_series_slice1_detailed_analysis.fig'];
%         saveas(gcf,savename);
%         results = [results BW_wo1_area 'area1 wash out1' BW_pers1_area 'area1 persistent1' BW_plat1_area 'area1 plateau1'];
%     elseif user_entry_slice == 2
%         savename = [save_folder '/time_series_slice2_detailed_analysis.png'];
%         saveas(gcf,savename)
%         savename = [save_folder '/time_series_slice2_detailed_analysis.fig'];
%         saveas(gcf,savename);
%         results = [results BW_wo2_area 'area2 wash out' BW_pers2_area 'area2 persistent' BW_plat2_area 'area2 plateau'];
%     elseif user_entry_slice == 3
%         savename = [save_folder '/time_series_slice3_detailed_analysis.png'];
%         saveas(gcf,savename)
%         savename = [save_folder '/time_series_slice3_detailed_analysis.fig'];
%         saveas(gcf,savename);
%         results = [results BW_wo3_area 'area3 wash out' BW_pers3_area 'area3 persistent' BW_plat3_area 'area3 plateau'];
%     end
    
    savename = [save_folder '/time_series_detailed_analysis.png'];
    saveas(gcf,savename)
    savename = [save_folder '/time_series_detailed_analysis.fig'];
    saveas(gcf,savename);
    BW_wo_volume = BW_wo_volume + BW_wo_area;
    BW_plat_volume = BW_plat_volume + BW_plat_area;
    BW_pers_volume = BW_pers_volume + BW_pers_area;
    if slice == nslices
        results = [results BW_wo_volume 'BW_wo_volume' BW_plat_volume 'BW_plat_volume' BW_pers_volume 'BW_pers_volume'];
    end

%     cut_off_rect_pos = [save_folder, '/cut_off_rect_pos.mat'];
%     load(cut_off_rect_pos, 'rect_pos');
%     x1 = round(rect_pos(1));
%     x2 = round(rect_pos(3));
%     y1 = round(rect_pos(2));
%     y2 = round(rect_pos(4));
%     my_axis = [x1 x1+x2 y1 y1+y2];

%     BW_wo_small = BW_wo((y1:y1+y2-1),x1:(x1+x2-1));
%     BW_pers_small = BW_pers((y1:y1+y2-1),x1:(x1+x2-1));
%     BW_plat_small = BW_plat((y1:y1+y2-1),x1:(x1+x2-1));
%     
%    if user_entry_slice == 1
%        BW_wo_small1 = BW_wo_small;
%        BW_pers_small1 = BW_pers_small;
%        BW_plat_small1 = BW_plat_small;
%    elseif user_entry_slice == 2
%        BW_wo_small2 = BW_wo_small;
%        BW_pers_small2 = BW_pers_small;
%        BW_plat_small2 = BW_plat_small;
%    elseif user_entry_slice == 3
%        BW_wo_small3 = BW_wo_small;
%        BW_pers_small3 = BW_pers_small;
%        BW_plat_small3 = BW_plat_small;
%    end


end

%     if flag1 ~= 1
%         BW_wo_small1 = zeros(512,512);
%         BW_pers_small1 = zeros(512,512);
%         BW_plat_small1 = zeros(512,512);
%     end
%     if flag2 ~= 1
%         BW_wo_small2 = zeros(512,512);
%         BW_pers_small2 = zeros(512,512);
%         BW_plat_small2 = zeros(512,512);
%     end
%     if flag3 ~= 1
%         BW_wo_small3 = zeros(512,512);
%         BW_pers_small3 = zeros(512,512);
%         BW_plat_small3 = zeros(512,512);
%     end
%     
%     % display of the spatial distribution of the three small areas
%     figure;
%     subplot(3,3,1); imshow(BW_wo_small1,[]); title(SID);
%     subplot(3,3,2); imshow(BW_pers_small1); 
%     subplot(3,3,3); imshow(BW_plat_small1);
%     subplot(3,3,4); imshow(BW_wo_small2,[]); 
%     subplot(3,3,5); imshow(BW_pers_small2); 
%     subplot(3,3,6); imshow(BW_plat_small2);
%     subplot(3,3,7); imshow(BW_wo_small3,[]); 
%     subplot(3,3,8); imshow(BW_pers_small3); 
%     subplot(3,3,9); imshow(BW_plat_small3);
% 
% 
%     savename = [save_folder '/detailed_areas_all.png'];
%     saveas(gcf,savename)
%     savename = [save_folder '/detailed_areas_all.fig'];
%     saveas(gcf,savename);


    % VOLUMES
lesion_vol = lesion_volume;
%     BW_wo_vol = BW_wo1_area + BW_wo2_area + BW_wo3_area;
    BW_wo_vol = BW_wo_vol/lesion_vol;
    BW_wo_vol_str = num2str(BW_wo_vol);
    wo_title = [SID '  ' 'wo ' BW_wo_vol_str];
%     BW_plat_vol = BW_plat1_area + BW_plat2_area + BW_plat3_area;
    BW_plat_vol = BW_plat_vol/lesion_vol;
    BW_plat_vol_str = num2str(BW_plat_vol);
    plat_title = ['plat ' BW_plat_vol_str];
%     BW_pers_vol = BW_pers1_area + BW_pers2_area + BW_pers3_area;
    BW_pers_vol = BW_pers_vol/lesion_vol;   
    BW_pers_vol_str = num2str(BW_pers_vol);
    pers_title = ['pers ' BW_pers_vol_str];
   
    

    %delta_vec = [0 delta1 delta1+delta2 delta1+delta2+delta3 delta1+delta2+delta3+delta4 delta1+delta2+delta3+delta4+delta5]/60;
    delta_vec = [1 2 3 4 5 6];
    figure;
    my_pixels2 = delta_vec;
%     FA_vec = [FA_wo_mean FA_plat_mean FA_pers_mean FA_lesion_mean_shifted FA_tissue_mean_shifted FA_tissue2_mean_shifted];
%     ADC_vec = [ADC_wo_mean ADC_plat_mean ADC_pers_mean ADC_lesion_mean_shifted ADC_tissue_mean_shifted ADC_tissue2_mean_shifted];
    wash_in_vec = [wash_in_wo_mean wash_in_plat_mean wash_in_pers_mean wash_in_ROI_DCE_mean wash_in_tissue_mean wash_in_tissue2_mean];
    wash_out_vec = [wash_out_wo_mean wash_out_plat_mean wash_out_pers_mean wash_out_ROI_DCE_mean wash_out_tissue_mean wash_out_tissue2_mean];
    subplot(4,3,1); plot(my_pixels2, DCEimage_wo_mean,'o'); title(wo_title);
    subplot(4,3,2); plot(my_pixels2, DCEimage_plat_mean,'+'); title(plat_title); 
    subplot(4,3,3); plot(my_pixels2, DCEimage_pers_mean,'x'); title(pers_title);
    subplot(4,3,4); plot(my_pixels2, DCEimage_lesion_mean,'^'); title('lesion');
    subplot(4,3,5); plot(my_pixels2, DCEimage_tissue_mean,'>'); title('tissue1');
    subplot(4,3,6); plot(my_pixels2, DCEimage_tissue2_mean,'<'); title('tissue2');
%     subplot(4,3,7); plot(my_pixels, FA_vec,'d'); title('FA (wo, plat, pers, l, t1, t2)');
%     subplot(4,3,8); plot(my_pixels, ADC_vec,'d');title('ADC (wo, plat, pers, l, t1, t2)');
    subplot(4,3,10); plot(my_pixels, wash_in_vec,'d');title('wash-in (wo, plat, pers, l, t1, t2)');
    subplot(4,3,11); plot(my_pixels, wash_out_vec,'d');title('wash-out (wo, plat, pers, l, t1, t2)');

    savename = [save_folder '/detailed_analysis_main.png'];
    saveas(gcf,savename)
    savename = [save_folder '/detailed_analysis_main.fig'];
    saveas(gcf,savename);
    


%% ADC & FA analysis with DCE ROI
% results = [results FA_pers_mean 'FA pers mean' FA_pers_std 'FA pers std' FA_plat_mean 'FA plat mean' FA_plat_std 'FA plat std' FA_wo_mean 'FA wo mean' FA_wo_std 'FA wo std' ];
% results = [results ADC_pers_mean 'ADC pers mean' ADC_pers_std 'ADC pers std' ADC_plat_mean 'ADC plat mean' ADC_plat_std 'ADC plat std' ADC_wo_mean 'ADC wo mean' ADC_wo_std 'ADC wo std'];

%% WASH-IN & WASH-OUT analysis with DCE ROI
results = [results wash_in_wo_mean 'wash-in wo mean' wash_in_wo_std 'wash-in wo std' wash_in_pers_mean 'wash-in pers mean' wash_in_pers_std 'wash-in pers std' wash_in_plat_mean 'wash-in plat mean' wash_in_plat_std 'wash-in plat std'];
results = [results wash_out_ROI_DCE_mean 'wash-out wo mean' wash_out_ROI_DCE_std 'wash-out wo std' wash_out_pers_mean 'wash-out pers mean' wash_out_pers_std 'wash-out pers std' wash_out_plat_mean 'wash-out plat mean' wash_out_plat_std 'wash-out plat std'];
results = [results wash_out_short_ROI_DCE_mean 'wash-out_short wo mean' wash_out_short_ROI_DCE_std 'wash-out_short wo std' wash_out_short_pers_mean 'wash-out_short pers mean' wash_out_short_pers_std 'wash-out_short pers std' wash_out_short_plat_mean 'wash-out_short plat mean' wash_out_short_plat_std 'wash-out_short plat std'];
results = [results DCEimage_wo_mean(1) 'DCE image no.1 wo mean' DCEimage_wo_mean(2) 'DCE image no.2 wo mean' DCEimage_wo_mean(3) 'DCE image no.3 wo mean' DCEimage_wo_mean(4) 'DCE image no.4 wo mean' DCEimage_wo_mean(5) 'DCE image no.5 wo mean' DCEimage_wo_mean(6) 'DCE image no.6 wo mean'];
results = [results DCEimage_pers_mean(1) 'DCE image no.1 pers mean' DCEimage_pers_mean(2) 'DCE image no.2 pers mean' DCEimage_pers_mean(3) 'DCE image no.3 pers mean' DCEimage_pers_mean(4) 'DCE image no.4 pers mean' DCEimage_pers_mean(5) 'DCE image no.5 pers mean' DCEimage_pers_mean(6) 'DCE image no.6 pers mean'];
results = [results DCEimage_plat_mean(1) 'DCE image no.1 plat mean' DCEimage_plat_mean(2) 'DCE image no.2 plat mean' DCEimage_plat_mean(3) 'DCE image no.3 plat mean' DCEimage_plat_mean(4) 'DCE image no.4 plat mean' DCEimage_plat_mean(5) 'DCE image no.5 plat mean' DCEimage_plat_mean(6) 'DCE image no.6 plat mean'];
results = [results DCEimage_wo_std(1) 'DCE image no.1 wo std' DCEimage_wo_std(2) 'DCE image no.2 wo std' DCEimage_wo_std(3) 'DCE image no.3 wo std' DCEimage_wo_std(4) 'DCE image no.4 wo std' DCEimage_wo_std(5) 'DCE image no.5 wo std' DCEimage_wo_std(6) 'DCE image no.6 wo std'];
results = [results DCEimage_pers_std(1) 'DCE image no.1 pers std' DCEimage_pers_std(2) 'DCE image no.2 pers std' DCEimage_pers_std(3) 'DCE image no.3 pers std' DCEimage_pers_std(4) 'DCE image no.4 pers std' DCEimage_pers_std(5) 'DCE image no.5 pers std' DCEimage_pers_std(6) 'DCE image no.6 pers std'];
results = [results DCEimage_plat_std(1) 'DCE image no.1 plat std' DCEimage_plat_std(2) 'DCE image no.2 plat std' DCEimage_plat_std(3) 'DCE image no.3 plat std' DCEimage_plat_std(4) 'DCE image no.4 plat std' DCEimage_plat_std(5) 'DCE image no.5 plat std' DCEimage_plat_std(6) 'DCE image no.6 plat std'];
%% END WASH-IN & WASH-OUT analysis with DCE ROI

results2 = cell(1);
for i =1:((size(results,2)-1)/2); results2(i,1)=results(1,2*i);end
for i =1:((size(results,2)-1)/2); results2(i,2)=results(1,1+2*i);end
results = results2;

savename = [save_folder '/detailed_analysis_results.mat'];
save(savename, 'results');
close all
savename = [save_folder '/wash_out_vector.mat'];
save(savename, 'wash_out_ROI_vec');

