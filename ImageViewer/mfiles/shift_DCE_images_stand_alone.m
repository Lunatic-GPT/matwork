function shift_DCE_images_stand_alone(params)

% SHIFT_DCE_IMAGES compensates the motion artifacts in the DCE series.
% Please provide the shifts with respect to the second DCE image of the DCE
% image series. The second image of the series is naturally not shifted!
% 'shift_vector' is a vector of the form [shift 2-1, shift 2-2 (= [0 0]), shift 2-3, shift 2-4,
% shift 2-5, shift 2-6], where each 'shift' consists of two values, i.e.
% shift in x-direction and shift in y-direction, according to the
% definitions given in the function 'shift_image'. The numbers in the
% vector given above refer to the number of the image in the DCE series.
% SID is the subject id, e.g. 'B49217'.

%load shift_gui

nslices = get(params,'Total number of slices');

%SID = get(params,'Enter the subject ID');

shift_vector(1) = get(params,'2-1x');
shift_vector(2) = get(params,'2-1y');
shift_vector(3) = 0;
shift_vector(4) = 0;
shift_vector(5) = get(params,'2-3x');
shift_vector(6) = get(params,'2-3y');
shift_vector(7) = get(params,'2-4x');
shift_vector(8) = get(params,'2-4y');
shift_vector(9) = get(params,'2-5x');
shift_vector(10) = get(params,'2-5y');
shift_vector(11) = get(params,'2-6x');
shift_vector(12) = get(params,'2-6y');

folder = get(params,'Select image directory name');
folder = [folder '/'];
%folder = ['/home/hahntobi/matlab/patient_images/' SID '/DCE_temp/'];
save_subfolder = folder;

total_slices = nslices;

for slices = 1:total_slices
    
    if slices == 1
        DCE_dir = uigetdir('/home/zong/breast','Please select the data folder');
        token = strtok(DCE_dir, 'X');
                user_entry_slice = 1;

    end
    user_entry_slice_str = num2str(user_entry_slice);
    DCE_dir = [token 'X' user_entry_slice_str];
    user_entry_slice = user_entry_slice + 1;
    %     final_BW_file = [DCE_dir, '/final_BWs.mat'];
    %     load(final_BW_file, 'inner_BW_large', 'tissue_BW_large', 'tissue2_BW_large');
    %     DCE_inner = inner_BW_large;
    %     DCE_inner_margin = bwmorph(DCE_inner, 'remove');
    %     DCE_tissue = tissue_BW_large;
    %     DCE_tissue2 = tissue2_BW_large;
    cut_off_rect_pos = [DCE_dir, '/cut_off_rect_pos.mat'];
    if exist(cut_off_rect_pos)
        load(cut_off_rect_pos, 'rect_pos');
        x1 = round(rect_pos(1));
        x2 = round(rect_pos(3));
        y1 = round(rect_pos(2));
        y2 = round(rect_pos(4));
    else
        x1 = 1;
        x2 = 511;
        y1 = 1;
        y2 = 511;
        display('I could not find the manually selected cut-off rectangle, therefore I am calculating the wash-out image in its full size. This may take longer than usual. ');
    end
    
    m = (slices-1)*6; % index of first image of series of specific slice
    slicenumber = num2str(slices);
    for i = 0:5

        shift_x = shift_vector(2*(i+1)-1);
        shift_y = shift_vector(2*(i+1));

        if (m+i+1) < 10
            index = num2str(m+i+1);
            index = ['0' index]; % to ensure the right ordering of the files
        else
            index = num2str(m+i+1);
        end
        image_name = [folder index];
        image = load(image_name, 'image');
        image = image.image;
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

        image = shift_image(shift_image(image,shift_y_abs,direc_y),shift_x_abs,direc_x);
        
        image_name_shifted = [image_name '_shifted'];
        save(image_name_shifted, 'image');
        save_filename = [image_name_shifted '.tiff'];
        imwrite(image/4000,save_filename,'tif');
        if i==0
            quot1 = image;
        elseif i==1 
            % WASH-IN image

            if exist('delta1')
                weight_factor = delta1;
            else weight_factor = 1;
            end

            quot2 = image;
            quot = zeros(size(image,1),size(image,2));
            quot(quot1~=0) = (quot2(quot1~=0)-quot1(quot1~=0))./(quot1(quot1~=0)*weight_factor);
            save_filename = [save_subfolder '/wash_in_slice' slicenumber];
            image = quot;
            wi = image;
            save(save_filename, 'image');
            save_filename = [save_filename '.tiff'];
            quot(quot>4)=4;
            quot = quot/max(max(quot)); % normalization of the wash-in image
            imwrite(quot,save_filename,'tif');
            


        end
    end

    %             % WASH-OUT
    for i=1:6
        DCE_number = 6*(slices - 1) + i;
        DCE_number = num2str(DCE_number);
        if (6*(slices-1)+i)<10
            DCE_number = ['0' DCE_number];
        end

            file = [save_subfolder DCE_number];

        load(file,'image');
        DCE(:,:,i) = image;
    end
    
    disp('Wash-out calculation --> This may take up to 3 minutes');
    image = zeros(512,512);
    my_array = cat(3,DCE(:,:,2), DCE(:,:,3), DCE(:,:,4), DCE(:,:,5), DCE(:,:,6));
    if exist('delta2')
        x_vec = [0 delta2 delta2+delta3 delta2+delta3+delta4 delta2+delta3+delta4+delta5]/3;
    else
        x_vec = [1 2 3 4 5];
    end

    for i = y1:(y1+y2)
        for j = x1:(x1+x2)
            temp = (squeeze(my_array(i,j,:)))';
            fit_vector=polyfit(x_vec,temp,1);
            slope = fit_vector(1);
            image(i,j) = slope;
        end
    end
    s_norm = (image/4000)/0.02;
    s_rad = atan(s_norm)*(2*90/pi);
    figure; imshow(s_rad,[])
    caxis([-90 90])
    jet_inv = jet;
    jet_inv = flipud(jet_inv);
    colormap(jet_inv)
    title('Wash-out');
    colorbar;
    savename = [save_subfolder '/wash_out' slicenumber];
    image = s_rad;
    wo = image;
    save(savename, 'image'); % save as .mat-file
    savename1 = [savename '.fig'];
    saveas(gcf,savename1)
    savename2 = [savename '_full_size.png'];
    saveas(gcf,savename2)
    image = image/max(max(image)); % normalization
    imwrite(image,savename2,'tif');
    close
    % END WASH-OUT
    
    
%     disp('Wash-out calculation --> This may take up to 3 minutes');
%     image = zeros(512,512);
%     my_array = cat(3,DCE(:,:,2), DCE(:,:,3), DCE(:,:,4), DCE(:,:,5), DCE(:,:,6));
%     if exist('delta2')
%         x_vec = [0 delta2 delta2+delta3 delta2+delta3+delta4 delta2+delta3+delta4+delta5]/3;
%     else
%         x_vec = [1 2 3 4 5];
%     end
% 
%     for i = 1: 512
%         for j = 1:512
%             temp = (squeeze(my_array(i,j,:)))';
%             fit_vector=polyfit(x_vec,temp,1);
%             slope = fit_vector(1);
%             image(i,j) = slope;
%         end
%     end
%     s_norm = (image/4000)/0.02;
%     s_rad = atan(s_norm)*(2*90/pi);
%     figure; imshow(s_rad,[])
%     caxis([-90 90])
%     jet_inv = jet;
%     jet_inv = flipud(jet_inv);
%     colormap(jet_inv)
%     title('Wash-out');
%     colorbar;
%     savename = [save_subfolder '/wash_out' slicenumber];
%     image = s_rad;
%     save(savename, 'image'); % save as .mat-file
%     savename1 = [savename '.fig'];
%     saveas(gcf,savename1)
%     savename2 = [savename '_full_size.png'];
%     saveas(gcf,savename2)
%     image = image/max(max(image)); % normalization
%     imwrite(image,savename2,'tif');
%     close
%     % END WASH-OUT

%     % WASH-OUT 'SHORT'
%     if exist('delta2')
%         weight_wash_out = (delta2+delta3+delta4+delta5)/3;
%     else
%         weight_wash_out = 4;
%     end
%     image = (DCE(:,:,6)-DCE(:,:,2))/weight_wash_out;
%     s_norm = (image/4000)/0.02;
%     s_rad = atan(s_norm)*(2*90/pi);
%     figure; imshow(s_rad,[])
%     caxis([-90 90])
%     jet_inv = jet;
%     jet_inv = flipud(jet_inv);
%     colormap(jet_inv)
%     title('Wash-out using only DCE image no.2 & no. 6');
%     colorbar;
%     savename = [save_subfolder '/wash_out_short' slicenumber];
%     image = s_rad;
%     save(savename, 'image'); % save as .mat-file
%     savename1 = [savename '.fig'];
%     saveas(gcf,savename1)
%     savename2 = [savename '_full_size.png'];
%     saveas(gcf,savename2)
%     image = image/max(max(image)); % normalization
%     imwrite(image,savename2,'tif');
%     close all
%     % END WASH-OUT 'SHORT'



end




