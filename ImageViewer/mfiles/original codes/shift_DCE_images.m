function shift_DCE_images(shift_vector,SID, number)

% SHIFT_DCE_IMAGES compensates the motion artifacts in the DCE series.
% Please provide the shifts with respect to the second DCE image of the DCE
% image series. The second image of the series is naturally not shifted!
% 'shift_vector' is a vector of the form [shift 2-1, shift 2-3, shift 2-4,
% shift 2-5, shift 2-6], where each 'shift' consists of two values, i.e.
% shift in x-direction and shift in y-direction, according to the
% definitions given in the function 'shift_image'. The numbers in the
% vector given above refer to the number of the image in the DCE series.
% SID is the subject id, e.g. 'B49217'.

folder = ['/home/hahntobi/matlab/patient_images/' SID '/DCE_temp/'];
save_subfolder = folder;

array_line = number;
load('time_points_array.mat');
nslices = time_points_array(array_line,8);
% total_slices = nslices;
total_slices = input('Please denote the total number of slices that have been downloaded: ');

disp('************* DCE SHIFT TOOL *******************');
disp('Please get the lesions first, based on the 2nd DCE image!');


for slices = 1:total_slices

    if slices == 1
        DCE_dir = uigetdir('/home/hahntobi/matlab/Data/','Please select the data folder');
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
    end

    %     my_axis = [x1 x1+x2 y1 y1+y2];

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
        
        image_name_shifted = [folder 'shifted_' index];
        save(image_name_shifted, 'image');
        save_filename = [image_name '.tiff'];
        imwrite(image/4000,save_filename,'tif')
        if i==0
            quot1 = image;
        elseif i==1 % save also the WASH-IN image

            if exist('delta1')
                weight_factor = delta1;
            else weight_factor = 1;
            end

            quot2 = image;
            quot = zeros(size(image,1),size(image,2));
            quot(quot1~=0) = (quot2(quot1~=0)-quot1(quot1~=0))./(quot1(quot1~=0)*weight_factor);
            save_filename = [save_subfolder '/wash_in_slice' slicenumber];
            image = quot;
            save(save_filename, 'image');
            save_filename = [save_filename '.tiff'];
            quot(quot>4)=4;
            quot = quot/max(max(quot)); % normalization of the wash-in image
            imwrite(quot,save_filename,'tif');
            %% Calculate Difference-images
            diff = (quot2 - quot1)./weight_factor;
            save_filename = [save_subfolder '/diff_slice' slicenumber];
            image = diff;
            save(save_filename, 'image');
            save_filename = [save_filename '.tiff'];
            figure; imshow(image,[]);
            saveas(gcf,save_filename);
            %% Calculate weighted Difference-images
            diff = (quot2 - quot1)./weight_factor;
            diff(diff<1)=1; %%% MENTION THAT !!!!! REASON THAT!!!!
            save_filename = [save_subfolder '/we_diff_slice' slicenumber];
            image = quot2.*diff;
            save(save_filename, 'image');
            save_filename = [save_filename '.tiff'];
            figure; imshow(quot2,[]);
            saveas(gcf,save_filename);
        end
    end

    %             % WASH-OUT
    DCE = zeros(512,512,6); % Adapt this in case the resolution or FOV changes
    for i=1:6
        DCE_number = 6*(slices - 1) + i;
        DCE_number = num2str(DCE_number);
        if (6*(slices-1)+i)<10
            DCE_number = ['0' DCE_number];
        end

            file = [save_subfolder 'shifted_' DCE_number];

        load(file,'image');
        DCE(:,:,i) = image;
    end
    
%     %% SPEED UP PART
%     h_rect = imrect(gca, [1 1 60 60]);
%     rect_api = iptgetapi(h_rect);
%     rect_api.setFixedAspectRatioMode(true)
%     disp('SPEED UP the wash-out calculation by selecting an area around the lesion. It MUST cover the lesion and the tissue areas!');
%     user_entry_rect_off = input('Do you want this to be the cropped image? (yes = 1, no = else)');
%     if user_entry_rect_off == 1
%         rect_pos = rect_api.getPosition();
%         rect_pos = round(rect_pos);
%     else
%         return;
%     end
    
%     for i=1:6
%         DCEimage = DCE(:,:,i);
%         tempimage = zeros(size(DCEimage,1),size(DCEimage,2));
%         tempimage(rect_pos(2):(rect_pos(2)+rect_pos(4)),rect_pos(1):(rect_pos(1)+rect_pos(3)))=DCEimage(rect_pos(2):(rect_pos(2)+rect_pos(4)),rect_pos(1):(rect_pos(1)+rect_pos(3)));
%         DCE(:,:,i)=tempimage;
%     end
    %% END SPEED UP
    
    
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
    save(savename, 'image'); % save as .mat-file
    savename1 = [savename '.fig'];
    saveas(gcf,savename1)
    savename2 = [savename '_full_size.png'];
    saveas(gcf,savename2)
    image = image/max(max(image)); % normalization
    imwrite(image,savename2,'tif');
    close
    % END WASH-OUT

    % WASH-OUT 'SHORT'
    if exist('delta2')
        weight_wash_out = (delta2+delta3+delta4+delta5)/3;
    else
        weight_wash_out = 4;
    end
    image = (DCE(:,:,6)-DCE(:,:,2))/weight_wash_out;
    s_norm = (image/4000)/0.02;
    s_rad = atan(s_norm)*(2*90/pi);
    figure; imshow(s_rad,[])
    caxis([-90 90])
    jet_inv = jet;
    jet_inv = flipud(jet_inv);
    colormap(jet_inv)
    title('Wash-out using only DCE image no.2 & no. 6');
    colorbar;
    savename = [save_subfolder '/wash_out_short' slicenumber];
    image = s_rad;
    save(savename, 'image'); % save as .mat-file
    savename1 = [savename '.fig'];
    saveas(gcf,savename1)
    savename2 = [savename '_full_size.png'];
    saveas(gcf,savename2)
    image = image/max(max(image)); % normalization
    imwrite(image,savename2,'tif');
    close all
    % END WASH-OUT 'SHORT'
end




