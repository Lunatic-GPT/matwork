function auto_shift_DCE_images_callback(params)

% SHIFT_DCE_IMAGES compensates the image shift due to motion in the DCE series.
% Please provide the shifts with respect to the second DCE image of the DCE
% image series. The second image of the series is naturally not shifted!
% 'shift_vector' is a vector of the form [shift 2-1, shift 2-2 (= [0 0]), shift 2-3, shift 2-4,
% shift 2-5, shift 2-6], where each 'shift' consists of two values, i.e.
% shift in x-direction and shift in y-direction, according to the
% definitions given in the function 'shift_image'. The numbers in the
% vector given above refer to the number of the image in the DCE series.

nslices = get(params,'Total number of slices');
fslice =get(params,'First slice');
folder = get(params,'Select image directory');  % folder for downloaded images. Also the folder for saving the shifted images.
mc =get(params,'Perform motion correction');
data_dir = get(params,'Select data directory');  % needed to get the crop region


images = cell(nslices,6);

if mc
  h_w = waitbar(0,'Motion correction');
end


for slice = 1:nslices 
   % read in the images. 
    for i = 1:6

        image_name = fullfile(folder,sprintf('%02d_%d',fslice+slice-1,i));
        img_temp = load(image_name, 'image');
        images{slice,i} = img_temp.image;
    end
    
    
    dir_struct = dir(fullfile(data_dir, ['*X',num2str(fslice+slice-1)]));
    
    if ~exist(fullfile(data_dir,dir_struct.name,'ROIs.mat'),'file')
        rect_pos = [1,1,size(images{slice,1},2)-1,size(images{slice,1},1)-1];
    else
        temp = load(fullfile(data_dir,dir_struct.name,'ROIs.mat'),'tissue2_BW_large_b');
        rect_pos = temp.tissue2_BW_large_b;
    end
    
    motionx = zeros(6,1);
    motiony = zeros(6,1);

    for i=1:6
        
        image_name = fullfile(folder,sprintf('%02d_%d',fslice+slice-1,i));
        image_name_shifted = [image_name '_shifted'];     
        tiff_filename = [image_name_shifted '.tiff'];
        
        
        if i==2
            image = images{slice,2};
            save(image_name_shifted, 'image');
            imwrite(image/4000,tiff_filename,'tif')
            continue;
        end
        
        if mc
            image1 = full2crop(images{slice,2},rect_pos);
            image2 =  full2crop(images{slice,i},rect_pos);
          [h,im_matched,theta,x,y]=image_registr_MI(uint16(image1),uint16(image2),0,10); 
          
          if i==1
            disp([mfilename,': Best alignment for slice ',num2str(fslice+slice-1),' (base t=2)']);
          end
          fprintf('t=%d: dx = %d, dy = %d.\n',i,x,y);
          images{slice,i} = circshift(images{slice,i},[y,x]);
          waitbar(((slice-1)*6+i)/nslices/6,h_w);
          motionx(i) = x;
          motiony(i) = y;
        end
        
         
        image = images{slice,i};
        save(image_name_shifted, 'image');
        imwrite(image/4000,tiff_filename,'tif');
        
    end
    
    if mc
        save(fullfile(data_dir,dir_struct.name,'motion_xy.mat'),'motionx','motiony');
    end
    
end
if mc
  delete(h_w);
end
%[h,im_matched, theta,x,y]=image_registr_MI(image1, image2, angle, step);


%% calculate wash in
       
  for slice = 1:nslices
            quot1 = images{slice,1};
            % WASH-IN image

            if exist('delta1','var')
                weight_factor = delta1;
            else weight_factor = 1;
            end

            quot2 = images{slice,2};
            sz = size(images{slice,2});
            quot = zeros(sz(1),sz(2));
            quot(quot1~=0) = (quot2(quot1~=0)-quot1(quot1~=0))./(quot1(quot1~=0)*weight_factor);
            image = quot;
            save_filename = fullfile(folder, sprintf('%02d_wash_in_norm_diff', fslice+slice-1));
            save(save_filename, 'image');
            quot(quot>4)=4;
            quot = quot/max(quot(:)); % normalization of the wash-in image
            imwrite(quot,[save_filename,'.tiff'],'tif');
            
            
            image = (images{slice,2}-images{slice,1})/90;  % changed from 80 to 90, 7-23-09, X.Z.
            image = atan(image)*(180/pi);
            save_filename = fullfile(folder, sprintf('%02d_wash_in', fslice+slice-1));
            save(save_filename, 'image');
            image = image/max(image(:));
            imwrite(image,[save_filename,'.tiff'],'tif');
            
            

   end
disp('Wash in calculation done');

%% calculate WASH-OUT
h = waitbar(0,'Wash out calculation');

for slice = 1:nslices
    
   
    % find the crop region.
    dir_struct = dir(fullfile(data_dir, ['*X',num2str(fslice+slice-1)]));
    
    cut_off_rect_pos = fullfile(data_dir,dir_struct.name, 'ROIs.mat');
    if exist(cut_off_rect_pos,'file')
        temp=load(cut_off_rect_pos, 'tissue2_BW_large_b');
        rect_pos = temp.tissue2_BW_large_b;
        x1 = round(rect_pos(1));
        x2 = round(rect_pos(3));
        y1 = round(rect_pos(2));
        y2 = round(rect_pos(4));
    else
        x1 = 1;
        x2 = 511;  % the right corner is x1+x2
        y1 = 1;
        y2 = 511; % the lower corner is y1+y2
        disp(['I could not find the manually selected crop region for slice ',num2str(fslice+slice-1)]);
        disp('Therefore I am calculating the wash-out image in its full size.');
        disp('This may take longer than usual. ');
    end
    

    sz = size(images{1,1});
    DCE = zeros(sz(1),sz(2),6);
    for i=1:6
        DCE(:,:,i) = images{slice,i};
    end
    
    disp(['Slice ',num2str(fslice+slice-1),': wash-out calculation...']);
    image = zeros(512,512);
    my_array = cat(3,DCE(:,:,2), DCE(:,:,3), DCE(:,:,4), DCE(:,:,5), DCE(:,:,6));
    if exist('delta2','var')
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
    s_norm = image/90;  % changed from 80 to 90, 7-23-09, X.Z.
    image = atan(s_norm)*(180/pi);
    
       
    savename = fullfile(folder, sprintf('%02d_wash_out', fslice+slice-1));
    save(savename, 'image'); % save as .mat-file
    
    image = image/max(max(image)); % normalization
    imwrite(image,[savename,'.tiff'],'tif');
    
    waitbar(slice/nslices,h);
     
end
 
delete(h);

disp('Done!');







