function detailed_DTI_analysis

% detailed_DTI_analysis uses an already created ROI and a certain cut-off
% criterion to analyze the ADC images in more detail, i.e. to calculate
% mean and std on the ADC images on different 'clusters'.

% open the specific folder
dir_name = uigetdir('/export/res/breast','Please select the patient data folder');
temp3 = [dir_name '/research_results_new_50_perc/slice3'];
temp2 = [dir_name '/research_results_new_50_perc/slice2'];
temp1 = [dir_name '/research_results_new_50_perc/slice1'];
flag1 = 0;
flag2 = 0;
flag3 = 0;

if exist(temp3)
    dir_data = temp3;
    flag3 = 1;
    temp = [dir_data '/final_BWs.mat'];
    load(temp);
    inner_BW_large3 = inner_BW_large;
elseif exist(temp2)
    dir_data = temp2;
elseif exist(temp1)
    dir_data = temp1;
end

if exist(temp2)
    flag2 = 1;
    temp = [temp2 '/final_BWs.mat'];
    load(temp);
    inner_BW_large2 = inner_BW_large;
end
if exist(temp1)
    flag1 = 1;
    temp = [temp1 '/final_BWs.mat'];
    load(temp);
    inner_BW_large1 = inner_BW_large;
end

figure
rec = gcf;

display('new patient');
mal = zeros(3,1);
ben = zeros(3,1);
plat= zeros(3,1);

shift_x = input('Please indicate the shift in the positive x-direction (Matlab convention): ');
shift_y = input('Please indicate the shift in the positive y-direction (Matlab convention): ');

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


col = -18.37; % cut-off value on the "left side" (standard in both cases 0)
cor = 25.67; % cut-off value on the "right side"

mal_mean = []; % '3D analysis'
ben_mean = [];
plat_mean = [];

if flag1 ==1
    temp = fullfile(dir_data, '/DCE_temp/wash_out1.mat');
    wash_out1 = load(temp, 'image');
    wash_out1 = wash_out1.image;
    wash_out1 = wash_out1.*inner_BW_large1;
    BW_wo1 = roicolor(wash_out1, -90,col).*inner_BW_large1;
    BW_plat1 = roicolor(wash_out1, col, cor).*inner_BW_large1;
    BW_pers1 = roicolor(wash_out1, cor, 90).*inner_BW_large1;
    BW_wo1_area = sum(sum(BW_wo1==1));
    BW_plat1_area = sum(sum(BW_plat1==1));
    BW_pers1_area = sum(sum(BW_pers1==1));
    
    %%% TEST
    temp = fullfile(dir_data, '/DCE_temp/wash_in_slice1.mat');
    wash_in1 = load(temp, 'image');
    wash_in1 = wash_in1.image;
    wash_in1 = wash_in1.*inner_BW_large1;
    temp_vec = wash_in1(:);
    temp_mean = mean(temp_vec(temp_vec~=0));
    wash_in1 = wash_in1./temp_mean;
    BW2_pers1 = roicolor(wash_in1, 0,0.9).*inner_BW_large1;
    BW2_plat1 = roicolor(wash_in1, 0.9, 1.0).*inner_BW_large1;
    BW2_wo1 = roicolor(wash_in1, 1.0, max(temp_vec)).*inner_BW_large1;
    BW2_wo1_area = sum(sum(BW2_wo1==1));
    BW2_plat1_area = sum(sum(BW2_plat1==1));
    BW2_pers1_area = sum(sum(BW2_pers1==1));
end

if flag2 ==1
    temp = fullfile(dir_data, '/DCE_temp/wash_out2.mat');
    wash_out2 = load(temp, 'image');
    wash_out2 = wash_out2.image;
    wash_out2 = wash_out2.*inner_BW_large2;
    BW_wo2 = roicolor(wash_out2, -90,col).*inner_BW_large2;
    BW_pers2 = roicolor(wash_out2, cor, 90).*inner_BW_large2;
    BW_plat2 = roicolor(wash_out2, col, cor).*inner_BW_large2;
    BW_wo2_area = sum(sum(BW_wo2==1));
    BW_plat2_area = sum(sum(BW_plat2==1));
    BW_pers2_area = sum(sum(BW_pers2==1));  
    
    %%% TEST
    temp = fullfile(dir_data, '/DCE_temp/wash_in_slice2.mat');
    wash_in2 = load(temp, 'image');
    wash_in2 = wash_in2.image;
    wash_in2 = wash_in2.*inner_BW_large2;
    temp_vec = wash_in2(:);
    temp_mean = mean(temp_vec(temp_vec~=0));
    wash_in2 = wash_in2./temp_mean;
    BW2_pers2 = roicolor(wash_in2, 0,0.9).*inner_BW_large2;
    BW2_plat2 = roicolor(wash_in2, 0.9, 1.0).*inner_BW_large2;
    BW2_wo2 = roicolor(wash_in2, 1.0, max(temp_vec)).*inner_BW_large2;
    BW2_wo2_area = sum(sum(BW2_wo2==1));
    BW2_plat2_area = sum(sum(BW2_plat2==1));
    BW2_pers2_area = sum(sum(BW2_pers2==1));
end

if flag3 ==1
    temp = fullfile(dir_data, '/DCE_temp/wash_out3.mat');
    wash_out3 = load(temp, 'image');
    wash_out3 = wash_out3.image;
    wash_out3 = wash_out3.*inner_BW_large3;
    BW_wo3 = roicolor(wash_out3, -90,col).*inner_BW_large3;
    BW_pers3 = roicolor(wash_out3, cor, 90).*inner_BW_large3;
    BW_plat3 = roicolor(wash_out3, col, cor).*inner_BW_large3;
    BW_wo3_area = sum(sum(BW_wo3==1));
    BW_plat3_area = sum(sum(BW_plat3==1));
    BW_pers3_area = sum(sum(BW_pers3==1));    
    
    %%% TEST
    temp = fullfile(dir_data, '/DCE_temp/wash_in_slice3.mat');
    wash_in3 = load(temp, 'image');
    wash_in3 = wash_in3.image;
    wash_in3 = wash_in3.*inner_BW_large3;
    temp_vec = wash_in3(:);
    temp_mean = mean(temp_vec(temp_vec~=0));
    wash_in3 = wash_in3./temp_mean;
    BW2_pers3 = roicolor(wash_in3, 0,0.9).*inner_BW_large3;
    BW2_plat3 = roicolor(wash_in3, 0.9, 1.0).*inner_BW_large3;
    BW2_wo3 = roicolor(wash_in3, 1.0, max(temp_vec)).*inner_BW_large3;
    BW2_wo3_area = sum(sum(BW2_wo3==1));
    BW2_plat3_area = sum(sum(BW2_plat3==1));
    BW2_pers3_area = sum(sum(BW2_pers3==1));
end

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
wash_in_ROI_vec = [];
wash_in_pers_vec = [];
wash_in_plat_vec = [];
wash_out_ROI_vec = [];
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
nslices = input('How many slices do you want to analyze? (up to 3)');

%%% SAVE FOLDER
save_folder_basic = dir_name;
% save_folder_basic = uigetdir('/export/res/breast','Please select the
% patient subject folder or the main save folder for the analysis results');
%save_folder2_basic = uigetdir('/home/hahntobi/matlab/results/', 'Please select a copy save folder outside of vali');
%%% END SAVE FOLDER

SID = input('Please enter the patient SID: ','s');

for slice = 1:nslices
    if slice == 1
        user_entry_slice = input('Assign the start slice number of the DCE series of interest: ');
    else
        user_entry_slice = user_entry_slice+1; % The case that only slice 1 and slice 3 are examined shouldn't happen!
    end
    user_entry_slice_str = num2str(user_entry_slice);
    slicen =num2str(slice);
    save_folder = [save_folder_basic '/research_results_new_50_perc/slice' user_entry_slice_str];
    %load_folder = [save_folder_basic '/research_results/slice' user_entry_slice_str];
    exist_test = exist(save_folder,'dir');
    if exist_test == 0
        mkdir(save_folder);
    end
    
    if user_entry_slice == 1
        BW_wo = BW_wo1;
        BW_pers = BW_pers1;
        BW_plat = BW_plat1;
        
        BW2_wo = BW2_wo1;
        BW2_pers = BW2_pers1;
        BW2_plat = BW2_plat1;
    elseif user_entry_slice == 2 
        BW_wo = BW_wo2;
        BW_pers = BW_pers2;
        BW_plat = BW_plat2;
        
        BW2_wo = BW2_wo2;
        BW2_pers = BW2_pers2;
        BW2_plat = BW2_plat2;
    elseif user_entry_slice == 3
        BW_wo = BW_wo3;
        BW_pers = BW_pers3;
        BW_plat = BW_plat3;
        
        BW2_wo = BW2_wo3;
        BW2_pers = BW2_pers3;
        BW2_plat = BW2_plat3;
    end
    
    DCE_image_dir = [dir_data '/DCE_temp'];
    
    %DCE_image_dir = ['/home/hahntobi/matlab/patient_images/' SID '/DCE_temp'];
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

    DCE_ref = double(DCE(:,:,1));
    DCE_ref(DCE_ref==0)=1; % setting zeros to one, we will forget these pixels later by using the 'no_zero_mask', just set so that the division by DCE_ref is possible
    for i=1:6
        DCE_temp = double(DCE(:,:,i)); % double-conversion to make sure that the calculation of the std is correct
        
%         no_zero_mask = ones(size(DCE_ref,1),size(DCE_ref,2));
%         no_zero_mask(DCE_ref==0)=(1-no_zero_mask(DCE_ref==0)); % is zero everywhere where DCE_ref == 0
%         no_zero_mask(DCE_temp==0)=0; % If DCE_TEMP is zero at one pixel, this pixel value is wrong!!!!


%         %% RELATIVE CHANGE CALULATION FOR TIME_SERIES
%         mask = BW_pers;
%         temp = no_zero_mask.*((DCE_temp-DCE_ref).*mask)./DCE_ref;
%         temp_vec = temp(mask~=0);
%         pers_avg(i)=mean(temp_vec);
%         pers_std(i)=std(temp_vec);
%         mask = BW_wo;
%         temp = no_zero_mask.*((DCE_temp-DCE_ref).*mask)./DCE_ref;
%         temp_vec = temp(mask~=0);
%         wo_avg(i)=mean(temp_vec);
%         wo_std(i)=std(temp_vec);
%         mask = BW_plat;
%         temp = no_zero_mask.*((DCE_temp-DCE_ref).*mask)./DCE_ref;
%         temp_vec = temp(mask~=0);
%         plat_avg(i)=mean(temp_vec);
%         plat_std(i)=std(temp_vec);
%         %% END RELATIVE CHANGE CALCULATION FOR TIME_SERIES

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
        
        
        %%% TEST
        mask = BW2_pers;
        temp = DCE_temp.*mask;
        temp_vec = temp(mask~=0);
        temp_vec(temp_vec==0)=[]; % regard the case that DCE_temp == 0
        r2_pers_avg(i)=mean(temp_vec);
        r2_pers_std(i)=std(temp_vec);
        mask = BW2_wo;
        temp = DCE_temp.*mask;
        temp_vec = temp(mask~=0);
        temp_vec(temp_vec==0)=[];
        r2_wo_avg(i)=mean(temp_vec);
        r2_wo_std(i)=std(temp_vec);
        mask = BW2_plat;
        temp = DCE_temp.*mask;
        temp_vec = temp(mask~=0);
        temp_vec(temp_vec==0)=[];
        r2_plat_avg(i)=mean(temp_vec);
        r2_plat_std(i)=std(temp_vec);
        %% END RAW VALUES CALCULATION FOR TIME_SERIES


        DCEimage_pers = DCE_temp.*BW_pers;
        DCEimage_pers_vec(i,(i-1)*512*512+1:i*512*512) = DCEimage_pers(:)';
        if slice == nslices
            temp = DCEimage_pers_vec(i,:);
            DCEimage_pers_mean(i) = mean(temp(temp~=0));
            DCEimage_pers_std(i) = std(temp(temp~=0));
        end

        DCEimage_plat = DCE_temp.*BW_plat;
        DCEimage_plat_vec(i,(i-1)*512*512+1:i*512*512) = DCEimage_plat(:)';
        if slice == nslices
            temp = DCEimage_plat_vec(i,:);
            DCEimage_plat_mean(i) = mean(temp(temp~=0));
            DCEimage_plat_std(i) = std(temp(temp~=0));
        end

        DCEimage_wo = DCE_temp.*BW_wo;
        DCEimage_wo_vec(i,(i-1)*512*512+1:i*512*512) = DCEimage_wo(:)';
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
    file = [dir_data '/ADC_temp/0' user_entry_slice_str '.mat'];
    %file = ['/home/hahntobi/matlab/patient_images/' SID '/ADC_temp/0' user_entry_slice_str '.mat'];
    load(file,'image');
    ADCimage = image;
    ADCimage = size_adapt(ADCimage);
    ADCimage = shift_image(shift_image(ADCimage,shift_y_abs,direc_y),shift_x_abs,direc_x);

    ADC_pers = ADCimage.*BW_pers;
    ADC_pers_vec = [ADC_pers_vec ADC_pers(:)];
    if slice == nslices
        ADC_pers_mean = mean(ADC_pers_vec(ADC_pers_vec~=0));
        ADC_pers_std = std(ADC_pers_vec(ADC_pers_vec~=0));
    end

    ADC_plat = ADCimage.*BW_plat;
    ADC_plat_vec = [ADC_plat_vec ADC_plat(:)];
    if slice == nslices
        ADC_plat_mean = mean(ADC_plat_vec(ADC_plat_vec~=0));
        ADC_plat_std = std(ADC_plat_vec(ADC_plat_vec~=0));
    end

    ADC_wo = ADCimage.*BW_wo;
    ADC_wo_vec = [ADC_wo_vec ADC_wo(:)];
    if slice == nslices
        ADC_wo_mean = mean(ADC_wo_vec(ADC_wo_vec~=0));
        ADC_wo_std = std(ADC_wo_vec(ADC_wo_vec~=0));
    end

    %% END ADC analysis with DCE ROI


    %%% FA image analysis start
    file = [dir_data '/FA_temp/0' user_entry_slice_str '.mat'];
    %file = ['/home/hahntobi/matlab/patient_images/' SID '/FA_temp/0' user_entry_slice_str '.mat'];
    load(file,'image');
    FAimage = image;
    FAimage = size_adapt(FAimage);
    FAimage = shift_image(shift_image(FAimage,shift_y_abs,direc_y),shift_x_abs,direc_x);
    
    FA_pers = FAimage.*BW_pers;
    FA_pers_vec = [FA_pers_vec FA_pers(:)];
    if slice == nslices
        FA_pers_mean = mean(FA_pers_vec(FA_pers_vec~=0));
        FA_pers_std = std(FA_pers_vec(FA_pers_vec~=0));
    end

    FA_plat = FAimage.*BW_plat;
    FA_plat_vec = [FA_plat_vec FA_plat(:)];
    if slice == nslices
        FA_plat_mean = mean(FA_plat_vec(FA_plat_vec~=0));
        FA_plat_std = std(FA_plat_vec(FA_plat_vec~=0));
    end

    FA_wo = FAimage.*BW_wo;
    FA_wo_vec = [FA_wo_vec FA_wo(:)];
    if slice == nslices
        FA_wo_mean = mean(FA_wo_vec(FA_wo_vec~=0));
        FA_wo_std = std(FA_wo_vec(FA_wo_vec~=0));
    end
    %% END FA image analysis

    %% WASH-IN Analysis
%     file = ['/home/hahntobi/matlab/patient_images/' SID '/DCE_temp/wash_in_slice' user_entry_slice_str '.mat'];
%     load(file,'image');
%     wash_in = image;
%     wash_in_ROI = wash_in.*BW_wo;
%     wash_in_ROI_single_vec = wash_in_ROI(:);
% 
%     wash_in_ROI_vec = [wash_in_ROI_vec wash_in_ROI_single_vec];
%     if slice == nslices
%         wash_in_ROI_DCE_mean = mean(wash_in_ROI_vec(wash_in_ROI_vec~=0));
%         wash_in_ROI_DCE_std = std(wash_in_ROI_vec(wash_in_ROI_vec~=0));
%     end
% 
%     wash_in_pers = wash_in.*BW_pers;
%     wash_in_pers_vec = [wash_in_pers_vec wash_in_pers(:)];
%     if slice == nslices
%         wash_in_pers_mean = mean(wash_in_pers_vec(wash_in_pers_vec~=0));
%         wash_in_pers_std = std(wash_in_pers_vec(wash_in_pers_vec~=0));
%     end
% 
%     wash_in_plat = wash_in.*BW_plat;
%     wash_in_plat_vec = [wash_in_plat_vec wash_in_plat(:)];
%     if slice == nslices
%         wash_in_plat_mean = mean(wash_in_plat_vec(wash_in_plat_vec~=0));
%         wash_in_plat_std = std(wash_in_plat_vec(wash_in_plat_vec~=0));
%     end
% 
%     %% END WASH-IN Analysis
% 
%     %% WASH-OUT Analysis
%     file = ['/home/hahntobi/matlab/patient_images/' SID '/DCE_temp/wash_out' user_entry_slice_str '.mat'];
%     load(file,'image');
%     wash_out = image;
%     wash_out_ROI = wash_out.*BW_wo;
%     wash_out_ROI_single_vec = wash_out_ROI(:);
%     
%     wash_out_ROI_vec = [wash_out_ROI_vec wash_out_ROI_single_vec];
%     if slice == nslices
%         wash_out_ROI_DCE_mean = mean(wash_out_ROI_vec(wash_out_ROI_vec~=0));
%         wash_out_ROI_DCE_std = std(wash_out_ROI_vec(wash_out_ROI_vec~=0));
%     end
% 
%     wash_out_pers = wash_out.*BW_pers;
%     wash_out_pers_vec = [wash_out_pers_vec wash_out_pers(:)];
%     if slice == nslices
%         wash_out_pers_mean = mean(wash_out_pers_vec(wash_out_pers_vec~=0));
%         wash_out_pers_std = std(wash_out_pers_vec(wash_out_pers_vec~=0));
%     end
% 
%     wash_out_plat = wash_out.*BW_plat;
%     wash_out_plat_vec = [wash_out_plat_vec wash_out_plat(:)];
%     if slice == nslices
%         wash_out_plat_mean = mean(wash_out_plat_vec(wash_out_plat_vec~=0));
%         wash_out_plat_std = std(wash_out_plat_vec(wash_out_plat_vec~=0));
%     end
% 
%     %% END WASH-OUT Analysis
% 
%     %% 'SHORT' WASH-OUT Analysis
%     file = ['/home/hahntobi/matlab/patient_images/' SID '/DCE_temp/wash_out_short' user_entry_slice_str '.mat'];
%     load(file,'image');
%     wash_out_short = image;
%     wash_out_short_ROI = wash_out_short.*BW_wo;
%     wash_out_short_ROI_single_vec = wash_out_short_ROI(:);
% 
%     wash_out_short_ROI_vec = [wash_out_short_ROI_vec wash_out_short_ROI_single_vec];
%     if slice == nslices
%         wash_out_short_ROI_DCE_mean = mean(wash_out_short_ROI_vec(wash_out_short_ROI_vec~=0));
%         wash_out_short_ROI_DCE_std = std(wash_out_short_ROI_vec(wash_out_short_ROI_vec~=0));
%     end
% 
%     wash_out_short_pers = wash_out_short.*BW_pers;
%     wash_out_short_pers_vec = [wash_out_short_pers_vec wash_out_short_pers(:)];
%     if slice == nslices
%         wash_out_short_pers_mean = mean(wash_out_short_pers_vec(wash_out_short_pers_vec~=0));
%         wash_out_short_pers_std = std(wash_out_short_pers_vec(wash_out_short_pers_vec~=0));
%     end
% 
%     wash_out_short_plat = wash_out_short.*BW_plat;
%     wash_out_short_plat_vec = [wash_out_short_plat_vec wash_out_short_plat(:)];
%     if slice == nslices
%         wash_out_short_plat_mean = mean(wash_out_short_plat_vec(wash_out_short_plat_vec~=0));
%         wash_out_short_plat_std = std(wash_out_short_plat_vec(wash_out_short_plat_vec~=0));
%     end
% 
%     %% END 'SHORT' WASH-OUT Analysis


    disp('____________________________________________________________________');
    disp('                      ANALYSIS RESULTS                              '); 

    % kinetic curves plots
    my_pixels = [1 2 3 4 5 6];
    figure;
    plot(my_pixels, r_pers_avg,'x', my_pixels, r_plat_avg,'+',my_pixels, r_wo_avg, 'o'); set(gca,'XTick',1:1:6); title('DCE series mean values for wash out, plateau and persistent areas');
    h = legend('pers average','plateau average','wo average','Location', 'Best'); set(h,'Interpreter','none')
    %%% (if-clause is not really necessary, it's automatically sorted into
    %%% the right slice-directory...)
    if user_entry_slice == 1
        savename = [save_folder '/time_series_slice1_detailed_analysis.png'];
        saveas(gcf,savename)
        savename = [save_folder '/time_series_slice1_detailed_analysis.fig'];
        saveas(gcf,savename);
        results = [results BW_wo1_area 'area1 wash out1' BW_pers1_area 'area1 persistent1' BW_plat1_area 'area1 plateau1'];
    elseif user_entry_slice == 2
        savename = [save_folder '/time_series_slice2_detailed_analysis.png'];
        saveas(gcf,savename)
        savename = [save_folder '/time_series_slice2_detailed_analysis.fig'];
        saveas(gcf,savename);
        results = [results BW_wo2_area 'area2 wash out' BW_pers2_area 'area2 persistent' BW_plat2_area 'area2 plateau'];
    elseif user_entry_slice == 3
        savename = [save_folder '/time_series_slice3_detailed_analysis.png'];
        saveas(gcf,savename)
        savename = [save_folder '/time_series_slice3_detailed_analysis.fig'];
        saveas(gcf,savename);
        results = [results BW_wo3_area 'area3 wash out' BW_pers3_area 'area3 persistent' BW_plat3_area 'area3 plateau'];
    end

    cut_off_rect_pos = [save_folder, '/cut_off_rect_pos.mat'];
    load(cut_off_rect_pos, 'rect_pos');
    x1 = round(rect_pos(1));
    x2 = round(rect_pos(3));
    y1 = round(rect_pos(2));
    y2 = round(rect_pos(4));
    my_axis = [x1 x1+x2 y1 y1+y2];

    BW_wo_small = BW_wo((y1:y1+y2-1),x1:(x1+x2-1));
    BW_pers_small = BW_pers((y1:y1+y2-1),x1:(x1+x2-1));
    BW_plat_small = BW_plat((y1:y1+y2-1),x1:(x1+x2-1));
    
    %%% TEST
    BW2_wo_small = BW2_wo((y1:y1+y2-1),x1:(x1+x2-1));
    BW2_pers_small = BW2_pers((y1:y1+y2-1),x1:(x1+x2-1));
    BW2_plat_small = BW2_plat((y1:y1+y2-1),x1:(x1+x2-1));

    % display of the different 'detailed analysis areas'
    figure;
    subplot(4,3,1); imshow(BW_wo_small,[]); 
    subplot(4,3,2); imshow(BW_pers_small); 
    subplot(4,3,3); imshow(BW_plat_small);
    subplot(4,3,4); plot(my_pixels, r_wo_avg,'o'); set(gca,'XTick',1:1:6); 
    subplot(4,3,5); plot(my_pixels, r_pers_avg,'x'); set(gca,'XTick',1:1:6); 
    subplot(4,3,6); plot(my_pixels, r_plat_avg,'+'); set(gca,'XTick',1:1:6); 
    subplot(4,3,7); imshow(BW2_wo_small,[]); 
    subplot(4,3,8); imshow(BW2_pers_small); 
    subplot(4,3,9); imshow(BW2_plat_small);
    subplot(4,3,10); plot(my_pixels, r2_wo_avg,'o'); set(gca,'XTick',1:1:6); 
    subplot(4,3,11); plot(my_pixels, r2_pers_avg,'x'); set(gca,'XTick',1:1:6); 
    subplot(4,3,12); plot(my_pixels, r2_plat_avg,'+'); set(gca,'XTick',1:1:6);     
    

    savename = [save_folder '/detailed_areas_wi_wo_comp.png'];
    saveas(gcf,savename)
    savename = [save_folder '/detailed_areas_wi_wo_comp.fig'];
    saveas(gcf,savename);
end

%% ADC & FA analysis with DCE ROI
results = [results FA_pers_mean 'FA pers mean' FA_pers_std 'FA pers std' FA_plat_mean 'FA plat mean' FA_plat_std 'FA plat std' FA_wo_mean 'FA wo mean' FA_wo_std 'FA wo std' ];
results = [results ADC_pers_mean 'ADC pers mean' ADC_pers_std 'ADC pers std' ADC_plat_mean 'ADC plat mean' ADC_plat_std 'ADC plat std' ADC_wo_mean 'ADC wo mean' ADC_wo_std 'ADC wo std'];

%% WASH-IN & WASH-OUT analysis with DCE ROI
%results = [results wash_in_ROI_DCE_mean 'wash-in wo mean' wash_in_ROI_DCE_std 'wash-in wo std' wash_in_pers_mean 'wash-in pers mean' wash_in_pers_std 'wash-in pers std' wash_in_plat_mean 'wash-in plat mean' wash_in_plat_std 'wash-in plat std'];
%results = [results wash_out_ROI_DCE_mean 'wash-out wo mean' wash_out_ROI_DCE_std 'wash-out wo std' wash_out_pers_mean 'wash-out pers mean' wash_out_pers_std 'wash-out pers std' wash_out_plat_mean 'wash-out plat mean' wash_out_plat_std 'wash-out plat std'];
%results = [results wash_out_short_ROI_DCE_mean 'wash-out_short wo mean' wash_out_short_ROI_DCE_std 'wash-out_short wo std' wash_out_short_pers_mean 'wash-out_short pers mean' wash_out_short_pers_std 'wash-out_short pers std' wash_out_short_plat_mean 'wash-out_short plat mean' wash_out_short_plat_std 'wash-out_short plat std'];
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

savename = [save_folder '/detailed_analysis_results_50perc.mat'];
save(savename, 'results');
close all

