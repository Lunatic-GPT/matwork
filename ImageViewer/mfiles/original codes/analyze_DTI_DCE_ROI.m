function analyze_DTI_DCE_ROI
% analyze_DTI_DCE_ROI lets you open DCE and DTI ROI Data files and analyzes
% the time behavior of the DCE image series inside the rim and core of the
% lesion and inside the tissue area around the DCE ROI.


load root_directory

user_entry = input('Do you want to start analyzing? (yes->1)');
if user_entry == 1    
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
    wi_lesion = [];
    wi_tissue = [];
    wi_tissue2 = [];
    
    results = cell(1); % result - cell
    nslices = input('How many slices do you want to analyze? (up to 3)');
    
    %%% SAVE FOLDER
    save_folder_basic = uigetdir('/export/res/breast','Please select the patient subject folder or the main save folder for the analysis results');
    %save_folder2_basic = uigetdir('/home/hahntobi/matlab/results/', 'Please select a copy save folder outside of vali');
    %%% END SAVE FOLDER
    
    %% IMAGE FOLDER (17/7/2008
    image_dir = uigetdir(root, 'Please select the image directory: (?\PatientID\DCEtemp)');
    %% END IMAGE FOLDER
    
    SID = input('Please enter the patient SID: ','s');
    savefoldername = input('Please denote the name of the folder where you want to save your date: ','s');
    
    for slice = 1:nslices
        if slice == 1
            user_entry_slice = input('Assign the number of the first slice: ');    
        else
            user_entry_slice = user_entry_slice+1; % The case that only slice 1 and slice 3 are examined shouldn't happen!
        end
        user_entry_slice_str = num2str(user_entry_slice);
        slicen =num2str(slice);        
        save_folder = [save_folder_basic '/' savefoldername '/slice' user_entry_slice_str];
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
            DCE_dir = uigetdir(root,'Please select the DCE data folder');
            token = strtok(DCE_dir, 'X');
        end
        DCE_dir = [token 'X' user_entry_slice_str];
        final_BW_file = [DCE_dir, '/final_BWs.mat'];
        load(final_BW_file, 'inner_BW_large', 'tissue_BW_large', 'tissue2_BW_large');
        DCE_inner = inner_BW_large;
        DCE_inner_margin = bwmorph(DCE_inner, 'remove');
        DCE_tissue = tissue_BW_large;
        DCE_tissue2 = tissue2_BW_large;
        
        %% REMOVAL OF NECROTIC AREAS, NEW! (7/19/2008)
        file = [image_dir '/' SID '/DCE_temp/wash_in_slice' user_entry_slice_str '.mat']; 
        load(file,'image');
        wash_in = image;  
        mean_tissue_uptake = mean(mean(wash_in.*DCE_tissue));
        std_tissue_uptake = std(std(wash_in.*DCE_tissue));
        uptake_threshold = mean_tissue_uptake + 2*std_tissue_uptake;
        uptake_mask = wash_in;
        uptake_mask(DCE_inner==0)=0;
        uptake_mask(uptake_mask<uptake_threshold)=0;
        uptake_mask(uptake_mask~=0)=1;
        uptake_mask = logical(uptake_mask);
        DCE_inner = DCE_inner.*uptake_mask;
        %% END REMOVAL OF NECROTIC AREAS
        

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

        DCE_image_dir = [image_dir '/' SID '/DCE_temp'];
        %DCE_image_dir = ['/home/hahntobi/matlab/patient_images_DTI/' SID '/DCE_temp'];
        for i=1:6 % six DCE images
            DCE_number = 6*(user_entry_slice - 1) + i;
            DCE_number = num2str(DCE_number);
            if (6*(user_entry_slice-1)+i)<10
                DCE_number = ['0' DCE_number];
            end
            file = [DCE_image_dir '/' DCE_number '_shifted'];
            load(file,'image');
            DCE(:,:,i) = image;
        end

        %% NEW WASH-IN (in degrees): created 6/28/2008
        image1 = DCE(:,:,1);
        image2 = DCE(:,:,2);
        new_image = image2-image1;
        new_image = new_image/80;
        new_image = atan(new_image)*(180/pi);
        temp = new_image.*inner_BW_large;
        temp2 = new_image.*tissue_BW_large;
        temp3 = new_image.*tissue2_BW_large;
        wi_lesion = [wi_lesion temp(:)];
        wi_tissue = [wi_tissue temp2(:)];
        wi_tissue2 = [wi_tissue2 temp3(:)];
        if slice == nslices
            wi_lesion_mean = mean(wi_lesion(wi_lesion~=0));
            wi_tissue_mean = mean(wi_tissue(wi_tissue~=0));
            wi_tissue2_mean = mean(wi_tissue2(wi_tissue2~=0));
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
        
        figure; imshow(DCE_long, 'InitialMagnification', 'fit'); colormap(jet); colorbar;  caxis([0 2000]);title('DCE series images with DCE tissue and lesion ROI');

        savename = [save_folder '/DCE_series.fig']; % save as fig...
        saveas(gcf,savename) 
        savename = [save_folder '/DCE_series.png']; % ... and png.
        saveas(gcf,savename) 
        
%         savename = [save_folder2 '/DCE_series.fig']; % save as fig...
%         saveas(gcf,savename) 
%         savename = [save_folder2 '/DCE_series.png']; % ... and png.
%         saveas(gcf,savename) 
%         
        figure; imshow(DCE_relative_change, 'InitialMagnification', 'fit'); colormap(jet); colorbar;  caxis([0 2]); title('DCE relative change images');
        
        savename = [save_folder '/DCE_rel_change_series.fig']; % save as fig...
        saveas(gcf,savename) 
        savename = [save_folder '/DCE_rel_change_series.png']; % ... and png.
        saveas(gcf,savename) 
        
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
            DCEimage_tissue_vec(i,(i-1)*512*512+1:i*512*512) = DCEimage_tissue(:)';
            if slice == nslices
                temp = DCEimage_tissue_vec(i,:);
                DCEimage_tissue_mean(i) = mean(temp(temp~=0));
                DCEimage_tissue_std(i) = std(temp(temp~=0));
            end

            DCEimage_tissue2 = DCE_temp.*DCE_tissue2;
            DCEimage_tissue2_vec(i,(i-1)*512*512+1:i*512*512) = DCEimage_tissue2(:)';
            if slice == nslices
                temp = DCEimage_tissue2_vec(i,:);
                DCEimage_tissue2_mean(i) = mean(temp(temp~=0));
                DCEimage_tissue2_std(i) = std(temp(temp~=0));
            end

            DCEimage_lesion = DCE_temp.*DCE_inner;
            DCEimage_lesion_vec(i,(i-1)*512*512+1:i*512*512) = DCEimage_lesion(:)';
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
        file = [image_dir '/' SID '/ADC_temp/0' user_entry_slice_str '.mat']; 
        %file = ['/home/hahntobi/matlab/patient_images_DTI/' SID '/ADC_temp/0' user_entry_slice_str '.mat']; 
        load(file,'image');
        ADCimage = image;
        ADCimage = size_adapt(ADCimage);
        
        ADC_temp = ADCimage((y1:y1+y2-1),x1:(x1+x2-1));
        %ADC_temp_DCE_ROI = ADC_temp.*(1-DCE_inner_margin_small); % should
        %be the same as below
        ADC_temp_DCE_ROI = ADC_temp.*sur_and_mar_inv;
        figure; imshow(ADC_temp_DCE_ROI, 'InitialMagnification', 'fit'); colorbar; caxis([0 3]); %colormap(jet)
        
        savename = [save_folder '/ADC_with_DCE_ROI.fig']; % save as fig...
        saveas(gcf,savename) 
        savename = [save_folder '/ADC_with_DCE_ROI.png']; % ... and png.
        saveas(gcf,savename) 
        
%         savename = [save_folder2 '/ADC_with_DCE_ROI.fig']; % save as fig...
%         saveas(gcf,savename) 
%         savename = [save_folder2 '/ADC_with_DCE_ROI.png']; % ... and png.
%         saveas(gcf,savename) 
            
        %% ADC analysis with DCE ROI
        ADC_ROI = ADCimage.*DCE_inner;
        ADC_ROI_vec = [ADC_ROI_vec ADC_ROI(:)];
        if slice == nslices % plot histogram after last slice
            figure; hist(ADC_ROI_vec(ADC_ROI_vec~=0),50);
            savename = [save_folder '/ADC_hist.png']; 
            saveas(gcf,savename)
            savename = [save_folder '/ADC_hist.fig']; 
            saveas(gcf,savename) 

%             savename = [save_folder2 '/ADC_hist.png']; 
%             saveas(gcf,savename)
%             savename = [save_folder2 '/ADC_hist.fig']; 
%             saveas(gcf,savename)
            %ADC_ROI_DCE_mean = mean(ADC_ROI_vec(ADC_ROI_vec~=0));
            %ADC_ROI_DCE_std = std(ADC_ROI_vec(ADC_ROI_vec~=0));
        end
        
        ADC_tissue = ADCimage.*DCE_tissue;
        ADC_tissue_vec = [ADC_tissue_vec ADC_tissue(:)];
        if slice == nslices
            ADC_tissue_mean = mean(ADC_tissue_vec(ADC_tissue_vec~=0));
            ADC_tissue_std = std(ADC_tissue_vec(ADC_tissue_vec~=0));
        end
        
        ADC_tissue2 = ADCimage.*DCE_tissue2;
        ADC_tissue2_vec = [ADC_tissue2_vec ADC_tissue2(:)];
        if slice == nslices
            ADC_tissue2_mean = mean(ADC_tissue2_vec(ADC_tissue2_vec~=0));
            ADC_tissue2_std = std(ADC_tissue2_vec(ADC_tissue2_vec~=0));
        end
        
        ADC_lesion = ADCimage.*DCE_inner;
        ADC_lesion_vec = [ADC_lesion_vec ADC_lesion(:)];
        if slice == nslices
            ADC_lesion_mean = mean(ADC_lesion_vec(ADC_lesion_vec~=0));
            ADC_lesion_std = std(ADC_lesion_vec(ADC_lesion_vec~=0));
        end

        %% END ADC analysis with DCE ROI
        
%             %% ADC analysis with DTI ROI
%             ADC_ROI = ADCimage.*DTI_inner;
%             ADC_ROI_vec = [ADC_ROI_vec ADC_ROI(:)];
%             if slice == 3 % plot histogram after third slice
%                 figure; hist(ADC_ROI_vec(ADC_ROI_vec~=0),50);
%                 savename = [save_folder '/ADC_hist.png']; 
%                 saveas(gcf,savename)
%                 savename = [save_folder '/ADC_hist.fig']; 
%                 saveas(gcf,savename) 
% 
%                 savename = [save_folder2 '/ADC_hist.png']; 
%                 saveas(gcf,savename)
%                 savename = [save_folder2 '/ADC_hist.fig']; 
%                 saveas(gcf,savename)
%             end
%             figure; subplot(2,2,1);imshow(ADC_ROI,[]);
%             axis(my_axis); title('DTI ROI assigned to the ADC image');
%             my_caxis = caxis;
%             temp = ADC_ROI(:);
%             temp = temp(temp~=0);
%             ADC_ROI_core=ADC_ROI;
%             ADC_ROI_core(ADC_ROI_core<mean(temp))=0;
%             ADC_ROI_rim = ADC_ROI - ADC_ROI_core;
%             ADC_mean_core = mean(mean(ADC_ROI_core(ADC_ROI_core~=0)));
%             ADC_mean_rim = mean(mean(ADC_ROI_rim(ADC_ROI_rim~=0)));
%             subplot(2,2,2);imshow(ADC_ROI_core,my_caxis);axis(my_axis); title('Core and outer part ADC image');
%             subplot(2,2,3);imshow(ADC_ROI_rim, my_caxis); axis(my_axis); title('Rim ADC image');
%             savename = [save_folder '/ADC_with_ROI.png']; 
%             saveas(gcf,savename) 
%             savename = [save_folder '/ADC_with_ROI.fig']; 
%             saveas(gcf,savename)     
% 
%             savename = [save_folder2 '/ADC_with_ROI.png']; 
%             saveas(gcf,savename) 
%             savename = [save_folder2 '/ADC_with_ROI.fig']; 
%             saveas(gcf,savename)
%             %% END ADC analysis with DTI ROI
            
        %%% ADC image analysis end
        
        %%% Difference % Weighted Difference image analysis start
%         file = [image_dir '/' SID '/DCE_temp/diff_slice' user_entry_slice_str '.mat']; 
%         %file = ['/home/hahntobi/matlab/patient_images_DTI/' SID '/DCE_temp/diff_slice' user_entry_slice_str '.mat']; 
%         load(file,'image');
%         diffimage = image;
%         diffimage_temp = diffimage((y1:y1+y2-1),x1:(x1+x2-1));
%         %diffimage_DCE_ROI = diffimage_temp.*sur_and_mar_inv; % WITH
%         %MARGINS
%         diffimage_DCE_ROI = diffimage_temp;
%         figure; imshow(diffimage_DCE_ROI, 'InitialMagnification', 'fit'); colorbar; colormap(jet);caxis([0 1600]);%colormap(jet);       
%         savename = [save_folder '/diff.fig']; % save as fig...
%         saveas(gcf,savename) 
%         savename = [save_folder '/diff.png']; % ... and png.
%         saveas(gcf,savename) 
%         
%         wdiffimage = diffimage.*DCE(:,:,2); % weighted Difference image
%         diffweight_temp = wdiffimage((y1:y1+y2-1),x1:(x1+x2-1));
%         %diffweight_DCE_ROI = diffweight_temp.*sur_and_mar_inv; % WITH
%         %MARGINS
%         diffweight_DCE_ROI = diffweight_temp;
%         figure; imshow(diffweight_DCE_ROI, 'InitialMagnification', 'fit'); colorbar; colormap(jet); caxis([0 2E6]);%colormap(jet);       
%         savename = [save_folder '/wdiff.fig']; % save as fig...
%         saveas(gcf,savename) 
%         savename = [save_folder '/wdiff.png']; % ... and png.
%         saveas(gcf,savename) 
%         
%         diff_tissue = diffimage.*DCE_tissue;
%         diff_tissue_vec = [diff_tissue_vec diff_tissue(:)];
%         if slice == nslices
%             diff_tissue_mean = mean(diff_tissue_vec(diff_tissue_vec~=0));
%             diff_tissue_std = std(diff_tissue_vec(diff_tissue_vec~=0));
%         end
%         
%         diff_tissue2 = diffimage.*DCE_tissue2;
%         diff_tissue2_vec = [diff_tissue2_vec diff_tissue2(:)];
%         if slice == nslices
%             diff_tissue2_mean = mean(diff_tissue2_vec(diff_tissue2_vec~=0));
%             diff_tissue2_std = std(diff_tissue2_vec(diff_tissue2_vec~=0));
%         end
%         
%         diff_lesion = diffimage.*DCE_inner;
%         diff_lesion_vec = [diff_lesion_vec diff_lesion(:)];
%         if slice == nslices
%             diff_lesion_mean = mean(diff_lesion_vec(diff_lesion_vec~=0));
%             diff_lesion_std = std(diff_lesion_vec(diff_lesion_vec~=0));
%         end
%         
%         wdiff_tissue = wdiffimage.*DCE_tissue;
%         wdiff_tissue_vec = [wdiff_tissue_vec wdiff_tissue(:)];
%         if slice == nslices
%             wdiff_tissue_mean = mean(wdiff_tissue_vec(wdiff_tissue_vec~=0));
%             wdiff_tissue_std = std(wdiff_tissue_vec(wdiff_tissue_vec~=0));
%         end
%         
%         wdiff_tissue2 = wdiffimage.*DCE_tissue2;
%         wdiff_tissue2_vec = [wdiff_tissue2_vec wdiff_tissue2(:)];
%         if slice == nslices
%             wdiff_tissue2_mean = mean(wdiff_tissue2_vec(wdiff_tissue2_vec~=0));
%             wdiff_tissue2_std = std(wdiff_tissue2_vec(wdiff_tissue2_vec~=0));
%         end
%         
%         wdiff_lesion = wdiffimage.*DCE_inner;
%         wdiff_lesion_vec = [wdiff_lesion_vec wdiff_lesion(:)];
%         if slice == nslices
%             wdiff_lesion_mean = mean(wdiff_lesion_vec(wdiff_lesion_vec~=0));
%             wdiff_lesion_std = std(wdiff_lesion_vec(wdiff_lesion_vec~=0));
%         end
%         %%% Difference % Weighted Difference image analysis end
%         
        
        %%% FA image analysis start
%         [filename,pathname] = uigetfile({'*.m;*.mat','MATLAB Files (*.m,*.mat)'},'Select the slice corresponding FA image','/home/hahntobi/matlab/patient_images_DTI/');
%         file = [pathname filename];
        file = [image_dir '/' SID '/FA_temp/0' user_entry_slice_str '.mat']; 
        %file = ['/home/hahntobi/matlab/patient_images_DTI/' SID '/FA_temp/0' user_entry_slice_str '.mat']; 
        load(file,'image');
        FAimage = image;
        FAimage = size_adapt(FAimage);
        
        FA_temp = FAimage((y1:y1+y2-1),x1:(x1+x2-1));
        %FA_temp_DCE_ROI = FA_temp.*(1-DCE_inner_margin_small);
        FA_temp_DCE_ROI = FA_temp.*sur_and_mar_inv;
        figure; imshow(FA_temp_DCE_ROI, 'InitialMagnification', 'fit'); colorbar; caxis([0 1.2]);%colormap(jet);
        
        savename = [save_folder '/FA_with_DCE_ROI.fig']; % save as fig...
        saveas(gcf,savename) 
        savename = [save_folder '/FA_with_DCE_ROI.png']; % ... and png.
        saveas(gcf,savename) 
%         
%         savename = [save_folder2 '/FA_with_DCE_ROI.fig']; % save as fig...
%         saveas(gcf,savename) 
%         savename = [save_folder2 '/FA_with_DCE_ROI.png']; % ... and png.
%         saveas(gcf,savename) 
        
        %% FA analysis with DCE ROI
        FA_ROI = FAimage.*DCE_inner;
        FA_ROI_vec = [FA_ROI_vec FA_ROI(:)];
        if slice == nslices % plot histogram after last slice
            figure; hist(FA_ROI_vec(FA_ROI_vec~=0),50);
            savename = [save_folder '/FA_hist.png']; 
            saveas(gcf,savename)
            savename = [save_folder '/FA_hist.fig']; 
            saveas(gcf,savename) 

%             savename = [save_folder2 '/FA_hist.png']; 
%             saveas(gcf,savename)
%             savename = [save_folder2 '/FA_hist.fig']; 
%             saveas(gcf,savename) 
            %FA_ROI_DCE_mean = mean(FA_ROI_vec(FA_ROI_vec~=0));
            %FA_ROI_DCE_std = std(FA_ROI_vec(FA_ROI_vec~=0));
        end
        
        FA_tissue = FAimage.*DCE_tissue;
        FA_tissue_vec = [FA_tissue_vec FA_tissue(:)];
        if slice == nslices
            FA_tissue_mean = mean(FA_tissue_vec(FA_tissue_vec~=0));
            FA_tissue_std = std(FA_tissue_vec(FA_tissue_vec~=0));
        end
        
        FA_tissue2 = FAimage.*DCE_tissue2;
        FA_tissue2_vec = [FA_tissue2_vec FA_tissue2(:)];
        if slice == nslices
            FA_tissue2_mean = mean(FA_tissue2_vec(FA_tissue2_vec~=0));
            FA_tissue2_std = std(FA_tissue2_vec(FA_tissue2_vec~=0));
        end
        
        FA_lesion = FAimage.*DCE_inner;
        FA_lesion_vec = [FA_lesion_vec FA_lesion(:)];
        if slice == nslices
            FA_lesion_mean = mean(FA_lesion_vec(FA_lesion_vec~=0));
            FA_lesion_std = std(FA_lesion_vec(FA_lesion_vec~=0));
        end
        %% END FA analysis with DCE ROI
        
        %% DTI analysis with DCE ROI
%         [filename,pathname] = uigetfile({'*.m;*.mat','MATLAB Files (*.m,*.mat)'},'Select the slice corresponding unweighted DTI image','/home/hahntobi/matlab/patient_images_DTI/');
%         file = [pathname filename];
        DTI_number = (user_entry_slice-1)*7+1;
        DTI_number_str = num2str(DTI_number);
        if user_entry_slice > 2
            file = [image_dir '/' SID '/DTI_temp/' DTI_number_str '.mat']; 
            %file = ['/home/hahntobi/matlab/patient_images_DTI/' SID '/DTI_temp/' DTI_number_str '.mat']; 
        else
            file = [image_dir '/' SID '/DTI_temp/0' DTI_number_str '.mat']; 
            %file = ['/home/hahntobi/matlab/patient_images_DTI/' SID '/DTI_temp/0' DTI_number_str '.mat'];
        end
        load(file,'image');
        unw_DTI = image; 
        unw_DTI = size_adapt(unw_DTI);
        unw_DTI_temp = unw_DTI((y1:y1+y2-1),x1:(x1+x2-1));
        unw_DTI_temp_DCE_ROI = unw_DTI_temp.*(1-DCE_inner_margin_small);
        figure; imshow(unw_DTI_temp_DCE_ROI, 'InitialMagnification', 'fit'); colorbar; caxis([0 800]); %colormap(jet);
        
        savename = [save_folder '/unw_DTI_with_DCE_ROI.fig']; % save as fig...
        saveas(gcf,savename) 
        savename = [save_folder '/unw_DTI_with_DCE_ROI.png']; % ... and png.
        saveas(gcf,savename) 
        
%         savename = [save_folder2 '/unw_DTI_with_DCE_ROI.fig']; % save as fig...
%         saveas(gcf,savename) 
%         savename = [save_folder2 '/unw_DTI_with_DCE_ROI.png']; % ... and png.
%         saveas(gcf,savename) 
%         
%         [filename,pathname] = uigetfile({'*.m;*.mat','MATLAB Files (*.m,*.mat)'},'Select the slice corresponding first weighted DTI image','/home/hahntobi/matlab/patient_images_DTI/');
%         file = [pathname filename];
        DTI_number = (user_entry_slice-1)*7+2;
        DTI_number_str = num2str(DTI_number);
%         if user_entry_slice > 2
%             file = ['/home/hahntobi/matlab/patient_images_DTI/' SID '/DTI_temp/' DTI_number_str '.mat']; 
%         else
%             file = ['/home/hahntobi/matlab/patient_images_DTI/' SID '/DTI_temp/0' DTI_number_str '.mat']; 
%         end
        if user_entry_slice > 2
            file = [image_dir '/' SID '/DTI_temp/' DTI_number_str '.mat']; 
            %file = ['/home/hahntobi/matlab/patient_images_DTI/' SID '/DTI_temp/' DTI_number_str '.mat']; 
        else
            file = [image_dir '/' SID '/DTI_temp/0' DTI_number_str '.mat']; 
            %file = ['/home/hahntobi/matlab/patient_images_DTI/' SID '/DTI_temp/0' DTI_number_str '.mat'];
        end
        load(file,'image');
        first_DTI = image;
        first_DTI = size_adapt(first_DTI);
        
        first_DTI_temp = first_DTI((y1:y1+y2-1),x1:(x1+x2-1));
        first_DTI_temp_DCE_ROI = first_DTI_temp.*(1-DCE_inner_margin_small);
        figure; imshow(first_DTI_temp_DCE_ROI, 'InitialMagnification', 'fit'); colorbar; caxis([0 400]);%colormap(jet);
        
        savename = [save_folder '/first_DTI_with_DCE_ROI.fig']; % save as fig...
        saveas(gcf,savename) 
        savename = [save_folder '/first_DTI_with_DCE_ROI.png']; % ... and png.
        saveas(gcf,savename) 
%         
%         savename = [save_folder2 '/first_DTI_with_DCE_ROI.fig']; % save as fig...
%         saveas(gcf,savename) 
%         savename = [save_folder2 '/first_DTI_with_DCE_ROI.png']; % ... and png.
%         saveas(gcf,savename) 
        
        DTI_unw_ROI = unw_DTI.*DCE_inner;
        DTI_unw_ROI_vec = [DTI_unw_ROI_vec DTI_unw_ROI(DCE_inner~=0)'];
        if slice == nslices
            DTI_unw_ROI_mean = mean(DTI_unw_ROI_vec);
            DTI_unw_ROI_std = std(DTI_unw_ROI_vec);            
        end

        DTI_first_ROI = first_DTI.*DCE_inner;
        DTI_first_ROI_vec = [DTI_first_ROI_vec DTI_first_ROI(:)];
        if slice == nslices
            DTI_first_ROI_mean = mean(DTI_first_ROI_vec(DTI_first_ROI_vec~=0));
            DTI_first_ROI_std = std(DTI_first_ROI_vec(DTI_first_ROI_vec~=0));            
        end


        file = [image_dir '/' SID '/DTI_temp/' 'slice_' user_entry_slice_str '_avg_DTI.mat'];
        %file = ['/home/hahntobi/matlab/patient_images_DTI/' SID '/DTI_temp/' 'slice_' user_entry_slice_str '_avg_DTI.mat'];
        load(file,'image');
        avg_DTI = image;
        avg_DTI = size_adapt(avg_DTI);
        
        avg_DTI_temp = avg_DTI((y1:y1+y2-1),x1:(x1+x2-1));
        avg_DTI_temp_DCE_ROI = avg_DTI_temp.*(1-DCE_inner_margin_small);
        figure; imshow(avg_DTI_temp_DCE_ROI, 'InitialMagnification', 'fit'); colorbar; caxis([0 400]);%colormap(jet);
        
        savename = [save_folder '/avg_DTI_with_DCE_ROI.fig']; % save as fig...
        saveas(gcf,savename) 
        savename = [save_folder '/avg_DTI_with_DCE_ROI.png']; % ... and png.
        saveas(gcf,savename) 
        
%         savename = [save_folder2 '/avg_DTI_with_DCE_ROI.fig']; % save as fig...
%         saveas(gcf,savename) 
%         savename = [save_folder2 '/avg_DTI_with_DCE_ROI.png']; % ... and png.
%         saveas(gcf,savename) 
        
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
%         [filename,pathname] = uigetfile({'*.m;*.mat','MATLAB Files (*.m,*.mat)'},'Select the slice corresponding WASH-IN image','/home/hahntobi/matlab/patient_images_DTI/');
%         file = [pathname filename];
        file = [image_dir '/' SID '/DCE_temp/wash_in_slice' user_entry_slice_str '.mat']; 
        %file = ['/home/hahntobi/matlab/patient_images_DTI/' SID '/DCE_temp/wash_in_slice' user_entry_slice_str '.mat']; 
        load(file,'image');
        wash_in = image;       
        wash_in_ROI = wash_in.*DCE_inner;
        wash_in_ROI_single_vec = wash_in_ROI(:);
        figure; hist(wash_in_ROI_single_vec(wash_in_ROI_single_vec~=0),50);
        savename = [save_folder '/wash_in_hist.png']; 
        saveas(gcf,savename)
        savename = [save_folder '/wash_in_hist.fig']; 
        saveas(gcf,savename) 
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
        figure; imshow(wash_in_temp, 'InitialMagnification', 'fit'); colormap(jet);colorbar; caxis([0 2])
        
        savename = [save_folder '/wash_in.fig']; % save as fig...
        saveas(gcf,savename) 
        savename = [save_folder '/wash_in.png']; % ... and png.
        saveas(gcf,savename) 
        
%         savename = [save_folder2 '/wash_in.fig']; % save as fig...
%         saveas(gcf,savename) 
%         savename = [save_folder2 '/wash_in.png']; % ... and png.
%         saveas(gcf,savename) 
        %% END WASH-IN Analysis
        
        %% WASH-OUT Analysis
        file = [image_dir '/' SID '/DCE_temp/wash_out' user_entry_slice_str '.mat']; 
        %file = ['/home/hahntobi/matlab/patient_images_DTI/' SID '/DCE_temp/wash_out' user_entry_slice_str '.mat']; 
        load(file,'image');
        wash_out = image;       
        wash_out_ROI = wash_out.*DCE_inner;
        wash_out_ROI_single_vec = wash_out_ROI(:);
        figure; hist(wash_out_ROI_single_vec(wash_out_ROI_single_vec~=0),50);
        savename = [save_folder '/wash_out_hist.png']; 
        saveas(gcf,savename)
        savename = [save_folder '/wash_out_hist.fig']; 
        saveas(gcf,savename) 
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
        figure; imshow(wash_out_temp, 'InitialMagnification', 'fit')
        caxis([-90 90])
        jet_inv = jet;
        jet_inv = flipud(jet_inv);
        colormap(jet_inv)
        colorbar
        
        savename = [save_folder '/wash_out.fig']; % save as fig...
        saveas(gcf,savename) 
        savename = [save_folder '/wash_out.png']; % ... and png.
        saveas(gcf,savename) 
        
%         savename = [save_folder2 '/wash_out.fig']; % save as fig...
%         saveas(gcf,savename) 
%         savename = [save_folder2 '/wash_out.png']; % ... and png.
%         saveas(gcf,savename) 
        %% END WASH-OUT Analysis
        
        %% 'SHORT' WASH-OUT Analysis
%         file = [image_dir '/' SID '/DCE_temp/wash_out_short' user_entry_slice_str '.mat']; 
%         %file = ['/home/hahntobi/matlab/patient_images_DTI/' SID '/DCE_temp/wash_out_short' user_entry_slice_str '.mat']; 
%         load(file,'image');
%         wash_out_short = image;       
%         wash_out_short_ROI = wash_out_short.*DCE_inner;
%         wash_out_short_ROI_single_vec = wash_out_short_ROI(:);
%         figure; hist(wash_out_short_ROI_single_vec(wash_out_short_ROI_single_vec~=0),50);
%         savename = [save_folder '/wash_out_short_hist.png']; 
%         saveas(gcf,savename)
%         savename = [save_folder '/wash_out_short_hist.fig']; 
%         saveas(gcf,savename) 
% %         savename = [save_folder2 '/wash_out_short_hist.png']; 
% %         saveas(gcf,savename)
% %         savename = [save_folder2 '/wash_out_short_hist.fig']; 
% %         saveas(gcf,savename) 
% 
%         wash_out_short_ROI_vec = [wash_out_short_ROI_vec wash_out_short_ROI_single_vec];
%         if slice == nslices
%             wash_out_short_ROI_DCE_mean = mean(wash_out_short_ROI_vec(wash_out_short_ROI_vec~=0));
%             wash_out_short_ROI_DCE_std = std(wash_out_short_ROI_vec(wash_out_short_ROI_vec~=0));        
%         end
% 
%         wash_out_short_tissue = wash_out_short.*DCE_tissue;
%         wash_out_short_tissue_vec = [wash_out_short_tissue_vec wash_out_short_tissue(:)];
%         if slice == nslices
%             wash_out_short_tissue_mean = mean(wash_out_short_tissue_vec(wash_out_short_tissue_vec~=0));
%             wash_out_short_tissue_std = std(wash_out_short_tissue_vec(wash_out_short_tissue_vec~=0));        
%         end
%         
%         wash_out_short_tissue2 = wash_out_short.*DCE_tissue2;
%         wash_out_short_tissue2_vec = [wash_out_short_tissue2_vec wash_out_short_tissue2(:)];
%         if slice == nslices
%             wash_out_short_tissue2_mean = mean(wash_out_short_tissue2_vec(wash_out_short_tissue2_vec~=0));
%             wash_out_short_tissue2_std = std(wash_out_short_tissue2_vec(wash_out_short_tissue2_vec~=0));        
%         end
%         
%         wash_out_short_temp  = wash_out_short((y1:y1+y2-1),x1:(x1+x2-1)); % = DCE_long(3*x2+1:4*x2,((2-1)*x2+1):((2-1)*x2+x2)) SHOULD BE THE SAME!!! 
%         figure; imshow(wash_out_short_temp, 'InitialMagnification', 'fit')
%         caxis([-90 90])
%         jet_inv = jet;
%         jet_inv = flipud(jet_inv);
%         colormap(jet_inv)
%         colorbar
%         
%         savename = [save_folder '/wash_out_short.fig']; % save as fig...
%         saveas(gcf,savename) 
%         savename = [save_folder '/wash_out_short.png']; % ... and png.
%         saveas(gcf,savename) 
        
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
        results = [results fin_area 'area of the DCE ROI'];
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
        figure; 
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

        plot(my_pixels, tissue_avg,'x', my_pixels, tissue2_avg,'+',my_pixels, lesion_avg, 'o'); set(gca,'XTick',1:1:6); title('peripheral areas and lesion average, relative change');
        h = legend('tissue average','tissue 2 average','lesion average','Location', 'Best'); set(h,'Interpreter','none')
        savename = [save_folder '/time_series.png']; 
        saveas(gcf,savename) 
        savename = [save_folder '/time_series.fig']; 
        saveas(gcf,savename)    
        
        plot(my_pixels, d_tissue_avg,'x', my_pixels, d_tissue2_avg,'+',my_pixels, d_lesion_avg, 'o'); set(gca,'XTick',1:1:6); title('peripheral areas and lesion average, difference image');
        h = legend('tissue average','tissue 2 average','lesion average','Location', 'Best'); set(h,'Interpreter','none')
        savename = [save_folder '/time_series_diff.png']; 
        saveas(gcf,savename) 
        savename = [save_folder '/time_series_diff.fig']; 
        saveas(gcf,savename)
        
        plot(my_pixels, wd_tissue_avg,'x', my_pixels, wd_tissue2_avg,'+',my_pixels, wd_lesion_avg, 'o'); set(gca,'XTick',1:1:6); title('peripheral areas and lesion average, weighted difference');
        h = legend('tissue average','tissue 2 average','lesion average','Location', 'Best'); set(h,'Interpreter','none')
        savename = [save_folder '/time_series_wdiff.png']; 
        saveas(gcf,savename) 
        savename = [save_folder '/time_series_wdiff.fig']; 
        saveas(gcf,savename)
        
        plot(my_pixels, r_tissue_avg,'x', my_pixels, r_tissue2_avg,'+',my_pixels, r_lesion_avg, 'o'); set(gca,'XTick',1:1:6); title('DCE series mean values');
        h = legend('tissue average','tissue 2 average','lesion average','Location', 'Best'); set(h,'Interpreter','none')
        savename = [save_folder '/time_series_basic.png']; 
        saveas(gcf,savename) 
        savename = [save_folder '/time_series_basic.fig']; 
        saveas(gcf,savename)
        
%         
%         savename = [save_folder2 '/time_series.png']; 
%         saveas(gcf,savename) 
%         savename = [save_folder2 '/time_series.fig']; 
%         saveas(gcf,savename)
        %my_time = clock;
        %name  = ['Research_results', num2str(my_time(2)),'_',num2str(my_time(3)),'_',num2str(my_time(1)),'_',num2str(my_time(4)),'_',num2str(my_time(5))];
        
        if slice ==1 
            %shift_question = input('Do you want to compensate a shift? (yes --> 1) ');
            shift_question=1;
        end
        if shift_question == 1
            %% ADC image analysis start
            shift_test = 0;
            figure;
            while shift_test == 0      
                if slice ==1 
                    shift_x = input('Please indicate the shift in the positive x-direction (Matlab convention): ');
                    shift_y = input('Please indicate the shift in the positive y-direction (Matlab convention): ');
                end
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
                ADC_temp = shift_image(shift_image(ADCimage,shift_y_abs,direc_y),shift_x_abs,direc_x);
                ADC_temp_small = ADC_temp((y1:y1+y2-1),x1:(x1+x2-1));
                FA_temp = shift_image(shift_image(FAimage,shift_y_abs,direc_y),shift_x_abs,direc_x);
                FA_temp_small = FA_temp((y1:y1+y2-1),x1:(x1+x2-1));
                unw_DTI_temp = shift_image(shift_image(unw_DTI,shift_y_abs,direc_y),shift_x_abs,direc_x);
                unw_DTI_temp = unw_DTI_temp((y1:y1+y2-1),x1:(x1+x2-1));
                first_DTI_temp = shift_image(shift_image(first_DTI,shift_y_abs,direc_y),shift_x_abs,direc_x);
                first_DTI_temp = first_DTI_temp((y1:y1+y2-1),x1:(x1+x2-1));
                avg_DTI_temp = shift_image(shift_image(avg_DTI,shift_y_abs,direc_y),shift_x_abs,direc_x);
                avg_DTI_temp = avg_DTI_temp((y1:y1+y2-1),x1:(x1+x2-1));
%                 ADC_temp_DCE_ROI = ADC_temp.*sur_and_mar_inv;
                subplot(5,2,1);imshow(ADC_temp_DCE_ROI.*sur_and_mar_inv, 'InitialMagnification', 'fit'); colorbar;  caxis([min(ADC_temp_DCE_ROI(:)) max(ADC_temp_DCE_ROI(:))]); %colormap(jet)     
                subplot(5,2,2);imshow(ADC_temp_small.*sur_and_mar_inv, 'InitialMagnification', 'fit'); colorbar; caxis([min(ADC_temp(:)) max(ADC_temp(:))]);
                subplot(5,2,3);imshow(FA_temp_DCE_ROI.*sur_and_mar_inv, 'InitialMagnification', 'fit'); colorbar; caxis([min(FA_temp_DCE_ROI(:)) max(FA_temp_DCE_ROI(:))]);
                subplot(5,2,4);imshow(FA_temp_small.*sur_and_mar_inv, 'InitialMagnification', 'fit'); colorbar; caxis([min(FA_temp(:)) max(FA_temp(:))]);
                subplot(5,2,5);imshow(unw_DTI_temp_DCE_ROI.*sur_and_mar_inv, 'InitialMagnification', 'fit'); colorbar; caxis([min(unw_DTI_temp_DCE_ROI(:)) max(unw_DTI_temp_DCE_ROI(:))]);
                subplot(5,2,6);imshow(unw_DTI_temp.*sur_and_mar_inv, 'InitialMagnification', 'fit'); colorbar; caxis([min(unw_DTI_temp(:)) max(unw_DTI_temp(:))]);
                subplot(5,2,7);imshow(first_DTI_temp_DCE_ROI.*sur_and_mar_inv, 'InitialMagnification', 'fit'); colorbar; caxis([min(first_DTI_temp_DCE_ROI(:)) max(first_DTI_temp_DCE_ROI(:))]);
                subplot(5,2,8);imshow(first_DTI_temp.*sur_and_mar_inv, 'InitialMagnification', 'fit'); colorbar; caxis([min(first_DTI_temp(:)) max(first_DTI_temp(:))]);
                subplot(5,2,9);imshow(avg_DTI_temp_DCE_ROI.*sur_and_mar_inv, 'InitialMagnification', 'fit'); colorbar; caxis([min(avg_DTI_temp_DCE_ROI(:)) max(avg_DTI_temp_DCE_ROI(:))]);
                subplot(5,2,10);imshow(avg_DTI_temp.*sur_and_mar_inv, 'InitialMagnification', 'fit'); colorbar; caxis([min(avg_DTI_temp(:)) max(avg_DTI_temp(:))]);
                if slice ==1
                    shift_test = input('Is this how you wanted the image to be shifted? (yes-->1)');
                else 
                    shift_test = 1; % the DTI/DCE should be the same for all slices
                end
            end
            
            shift_word1 = num2str(shift_x_abs);
            if shift_x < 0
                shift_word1 = ['minus_' shift_word1];
            end

            shift_word2 = num2str(shift_y_abs);
            if shift_y < 0
                shift_word2 = ['minus_' shift_word2];
            end
            
            savename = [save_folder '/DTI_shift_x_' shift_word1 '_y_' shift_word2 '.fig']; % save as fig...
            saveas(gcf,savename)
            savename = [save_folder '/DTI_shift_x_' shift_word1 '_y_' shift_word2 '.png']; % ... and png.
            saveas(gcf,savename)

%             savename = [save_folder2 '/DTI_shift_x_' shift_word1 '_y_' shift_word2 '.fig']; % save as fig...
%             saveas(gcf,savename)
%             savename = [save_folder2 '/DTI_shift_x_' shift_word1 '_y_' shift_word2 '.png']; % ... and png.
%             saveas(gcf,savename)
%             
            figure; imshow(ADC_temp_small,'InitialMagnification', 'fit'); colorbar; caxis([0 3]);
            
            savename = [save_folder '/ADC_with_DCE_ROI_shifted.fig']; % save as fig...
            saveas(gcf,savename)
            savename = [save_folder '/ADC_with_DCE_ROI_shifted.png']; % ... and png.
            saveas(gcf,savename)
% 
%             savename = [save_folder2 '/ADC_with_DCE_ROI_shifted.fig']; % save as fig...
%             saveas(gcf,savename)
%             savename = [save_folder2 '/ADC_with_DCE_ROI_shifted.png']; % ... and png.
%             saveas(gcf,savename)

            %% ADC analysis with DCE ROI


            ADC_tissue_shifted = ADC_temp.*DCE_tissue;
            ADC_tissue_vec_shifted = [ADC_tissue_vec_shifted ADC_tissue_shifted(:)];
            if slice == nslices
                ADC_tissue_mean_shifted = mean(ADC_tissue_vec_shifted(ADC_tissue_vec_shifted~=0));
                ADC_tissue_std_shifted = std(ADC_tissue_vec_shifted(ADC_tissue_vec_shifted~=0));
            end
            
            ADC_tissue2_shifted = ADC_temp.*DCE_tissue2;
            ADC_tissue2_vec_shifted = [ADC_tissue2_vec_shifted ADC_tissue2_shifted(:)];
            if slice == nslices
                ADC_tissue2_mean_shifted = mean(ADC_tissue2_vec_shifted(ADC_tissue2_vec_shifted~=0));
                ADC_tissue2_std_shifted = std(ADC_tissue2_vec_shifted(ADC_tissue2_vec_shifted~=0));
            end
            
            ADC_lesion_shifted = ADC_temp.*DCE_inner;
            ADC_lesion_vec_shifted = [ADC_lesion_vec_shifted ADC_lesion_shifted(:)];
            if slice == nslices
                ADC_lesion_mean_shifted = mean(ADC_lesion_vec_shifted(ADC_lesion_vec_shifted~=0));
                ADC_lesion_std_shifted = std(ADC_lesion_vec_shifted(ADC_lesion_vec_shifted~=0));
            end

            %% END ADC analysis with DCE ROI

            %%% FA image analysis start
            %         [filename,pathname] = uigetfile({'*.m;*.mat','MATLAB Files (*.m,*.mat)'},'Select the slice corresponding FA image','/home/hahntobi/matlab/patient_images_DTI/');
            %         file = [pathname filename];
            

            figure; imshow(FA_temp_small.*sur_and_mar_inv, 'InitialMagnification', 'fit'); colorbar; caxis([0 1.2]); %colormap(jet)

            savename = [save_folder '/FA_with_DCE_ROI_shifted.fig']; % save as fig...
            saveas(gcf,savename)
            savename = [save_folder '/FA_with_DCE_ROI_shifted.png']; % ... and png.
            saveas(gcf,savename)

%             savename = [save_folder2 '/FA_with_DCE_ROI_shifted.fig']; % save as fig...
%             saveas(gcf,savename)
%             savename = [save_folder2 '/FA_with_DCE_ROI_shifted.png']; % ... and png.
%             saveas(gcf,savename)

            FA_tissue_shifted = FA_temp.*DCE_tissue;
            FA_tissue_vec_shifted = [FA_tissue_vec_shifted FA_tissue_shifted(:)];
            if slice == nslices
                FA_tissue_mean_shifted = mean(FA_tissue_vec_shifted(FA_tissue_vec_shifted~=0));
                FA_tissue_std_shifted = std(FA_tissue_vec_shifted(FA_tissue_vec_shifted~=0));
            end
            
            FA_tissue2_shifted = FA_temp.*DCE_tissue2;
            FA_tissue2_vec_shifted = [FA_tissue2_vec_shifted FA_tissue2_shifted(:)];
            if slice == nslices
                FA_tissue2_mean_shifted = mean(FA_tissue2_vec_shifted(FA_tissue2_vec_shifted~=0));
                FA_tissue2_std_shifted = std(FA_tissue2_vec_shifted(FA_tissue2_vec_shifted~=0));
            end

            FA_lesion_shifted = FA_temp.*DCE_inner;
            FA_lesion_vec_shifted = [FA_lesion_vec_shifted FA_lesion_shifted(:)];
            if slice == nslices
                FA_lesion_mean_shifted = mean(FA_lesion_vec_shifted(FA_lesion_vec_shifted~=0));
                FA_lesion_std_shifted = std(FA_lesion_vec_shifted(FA_lesion_vec_shifted~=0));
            end
        end

    end
        %% ADC & FA analysis with DCE ROI
            %results = [results FA_ROI_DCE_mean 'FA_ROI_DCE_mean' FA_ROI_DCE_std 'FA_ROI_DCE_std' ADC_ROI_DCE_mean 'ADC_ROI_DCE_mean' ADC_ROI_DCE_std 'ADC_ROI_DCE_std'];
            %FA_ROI_DCE_mean_str = num2str(FA_ROI_DCE_mean);
            %FA_ROI_DCE_std_str = num2str(FA_ROI_DCE_std);
            %ADC_ROI_DCE_mean_str = num2str(ADC_ROI_DCE_mean);
            %ADC_ROI_DCE_std_str = num2str(ADC_ROI_DCE_std);
            %message = ['FA image with DCE ROI: ' FA_ROI_DCE_mean_str ' +- ' FA_ROI_DCE_std_str '. ADC image with DCE ROI: ' ADC_ROI_DCE_mean_str ' +- ' ADC_ROI_DCE_std_str];
            %disp(message);  
        %results = [results FA_tissue_mean 'FA tissue mean' FA_tissue_std 'FA tissue std' FA_tissue2_mean 'FA tissue 2 mean' FA_tissue2_std 'FA tissue 2 std' FA_lesion_mean 'FA lesion mean' FA_lesion_std 'FA lesion std' ];
        %results = [results ADC_tissue_mean 'ADC tissue mean' ADC_tissue_std 'ADC tissue std' ADC_tissue2_mean 'ADC tissue 2 mean' ADC_tissue2_std 'ADC tissue 2 std' ADC_lesion_mean 'ADC lesion mean' ADC_lesion_std 'ADC lesion std'];
        %% END ADC & FA analysis with DCE ROI
        
        %% DTI analysis with DCE ROI
        %results = [results DTI_unw_ROI_mean 'DTI_unweighted_DCE_mean' DTI_first_ROI_mean 'DTI_first_DCE_mean'];
        %% END DTI analysis with DCE ROI
        
        %% WASH-IN & WASH-OUT analysis with DCE ROI
        results = [results wash_in_ROI_DCE_mean 'wash-in lesion mean' wash_in_ROI_DCE_std 'wash-in lesion std' wash_in_tissue_mean 'wash-in tissue mean' wash_in_tissue_std 'wash-in tissue std' wash_in_tissue2_mean 'wash-in tissue 2 mean' wash_in_tissue2_std 'wash-in tissue 2 std'];
        results = [results wash_out_ROI_DCE_mean 'wash-out lesion mean' wash_out_ROI_DCE_std 'wash-out lesion std' wash_out_tissue_mean 'wash-out tissue mean' wash_out_tissue_std 'wash-out tissue std' wash_out_tissue2_mean 'wash-out tissue 2 mean' wash_out_tissue2_std 'wash-out tissue 2 std'];
        %results = [results wash_out_short_ROI_DCE_mean 'wash-out_short lesion mean' wash_out_short_ROI_DCE_std 'wash-out_short lesion std' wash_out_short_tissue_mean 'wash-out_short tissue mean' wash_out_short_tissue_std 'wash-out_short tissue std' wash_out_short_tissue2_mean 'wash-out_short tissue 2 mean' wash_out_short_tissue2_std 'wash-out_short tissue 2 std'];
        results = [results DCEimage_lesion_mean(1) 'DCE image no.1 lesion mean' DCEimage_lesion_mean(2) 'DCE image no.2 lesion mean' DCEimage_lesion_mean(3) 'DCE image no.3 lesion mean' DCEimage_lesion_mean(4) 'DCE image no.4 lesion mean' DCEimage_lesion_mean(5) 'DCE image no.5 lesion mean' DCEimage_lesion_mean(6) 'DCE image no.6 lesion mean'];
        results = [results DCEimage_tissue_mean(1) 'DCE image no.1 tissue mean' DCEimage_tissue_mean(2) 'DCE image no.2 tissue mean' DCEimage_tissue_mean(3) 'DCE image no.3 tissue mean' DCEimage_tissue_mean(4) 'DCE image no.4 tissue mean' DCEimage_tissue_mean(5) 'DCE image no.5 tissue mean' DCEimage_tissue_mean(6) 'DCE image no.6 tissue mean'];
        results = [results DCEimage_tissue2_mean(1) 'DCE image no.1 tissue2 mean' DCEimage_tissue2_mean(2) 'DCE image no.2 tissue2 mean' DCEimage_tissue2_mean(3) 'DCE image no.3 tissue2 mean' DCEimage_tissue2_mean(4) 'DCE image no.4 tissue2 mean' DCEimage_tissue2_mean(5) 'DCE image no.5 tissue2 mean' DCEimage_tissue2_mean(6) 'DCE image no.6 tissue2 mean'];
        results = [results DCEimage_lesion_std(1) 'DCE image no.1 lesion std' DCEimage_lesion_std(2) 'DCE image no.2 lesion std' DCEimage_lesion_std(3) 'DCE image no.3 lesion std' DCEimage_lesion_std(4) 'DCE image no.4 lesion std' DCEimage_lesion_std(5) 'DCE image no.5 lesion std' DCEimage_lesion_std(6) 'DCE image no.6 lesion std'];
        results = [results DCEimage_tissue_std(1) 'DCE image no.1 tissue std' DCEimage_tissue_std(2) 'DCE image no.2 tissue std' DCEimage_tissue_std(3) 'DCE image no.3 tissue std' DCEimage_tissue_std(4) 'DCE image no.4 tissue std' DCEimage_tissue_std(5) 'DCE image no.5 tissue std' DCEimage_tissue_std(6) 'DCE image no.6 tissue std'];
        results = [results DCEimage_tissue2_std(1) 'DCE image no.1 tissue2 std' DCEimage_tissue2_std(2) 'DCE image no.2 tissue2 std' DCEimage_tissue2_std(3) 'DCE image no.3 tissue2 std' DCEimage_tissue2_std(4) 'DCE image no.4 tissue2 std' DCEimage_tissue2_std(5) 'DCE image no.5 tissue2 std' DCEimage_tissue2_std(6) 'DCE image no.6 tissue2 std'];
        %results = [results diff_lesion_mean 'Difference image lesion mean' diff_tissue_mean 'Difference image tissue mean' diff_tissue2_mean 'Difference image tissue2 mean'];
        %results = [results diff_lesion_std 'Difference image lesion std' diff_tissue_std 'Difference image tissue std' diff_tissue2_std 'Difference image tissue2 std'];
        %results = [results wdiff_lesion_mean 'Weighted difference image lesion mean' wdiff_tissue_mean 'Weighted difference image tissue mean' wdiff_tissue2_mean 'Weighted difference image tissue2 mean'];
        %results = [results wdiff_lesion_std 'Weighted difference image lesion std' wdiff_tissue_std 'Weighted difference image tissue std' wdiff_tissue2_std 'Weighted difference image tissue2 std'];
        %% END WASH-IN & WASH-OUT analysis with DCE ROI
        
        
        %% append the optional shifted ADC/FA values
        if shift_question == 1
            results = [results FA_tissue_mean_shifted 'shifted FA tissue mean' FA_tissue_std_shifted 'shifted FA tissue std' FA_tissue2_mean_shifted 'shifted FA tissue 2 mean' FA_tissue2_std_shifted 'shifted FA tissue 2 std' FA_lesion_mean_shifted 'shifted FA lesion mean' FA_lesion_std_shifted 'shifted FA lesion std' ];
            results = [results ADC_tissue_mean_shifted 'shifted ADC tissue mean' ADC_tissue_std_shifted 'shifted ADC tissue std' ADC_tissue2_mean_shifted 'shifted ADC tissue 2 mean' ADC_tissue2_std_shifted 'shifted ADC tissue 2 std' ADC_lesion_mean_shifted 'shifted ADC lesion mean' ADC_lesion_std_shifted 'shifted ADC lesion std' ];
        end
        %% end appended shifted ADC/FA values
        
        %% append new wash-in values
            results = [results wi_lesion_mean 'wi_lesion_mean' wi_tissue_mean 'wi_tissue_mean' wi_tissue2_mean 'wi_tissue2_mean'];
        %% end append new wash-in values
        results2 = cell(1);
        for i =1:((size(results,2)-1)/2); results2(i,1)=results(1,2*i);end
        for i =1:((size(results,2)-1)/2); results2(i,2)=results(1,1+2*i);end
        results = results2;
        
        savename = [save_folder '/analysis_results.mat']; 
        save(savename, 'results');
%         savename = [save_folder2 '/analysis_results.mat'];         
%         save(savename, 'results');
        copy_folder = [image_dir '/' SID '/'];
        %copy_folder = ['/home/hahntobi/matlab/patient_images_DTI/' SID '/'];
        copyfile(copy_folder,save_folder); % copy all images
%         copyfile(copy_folder,save_folder2);
        
    close all
end