function compare_DCE_DTI_poss_shift
% Compare_DCE_DTI_poss_shift compares DCE and DTI images and determines if
% there is a possible (overall) shift between these images.

close all

DTI = zeros(256,256,7);

load('~hahntobi/matlab/work_images/DTI_temp/01', 'image');
DTI(:,:,1) = image;
load('~hahntobi/matlab/work_images/DTI_temp/02', 'image');
DTI(:,:,2) = image;
load('~hahntobi/matlab/work_images/DTI_temp/03', 'image');
DTI(:,:,3) = image;
load('~hahntobi/matlab/work_images/DTI_temp/04', 'image');
DTI(:,:,4) = image;
load('~hahntobi/matlab/work_images/DTI_temp/05', 'image');
DTI(:,:,5) = image;
load('~hahntobi/matlab/work_images/DTI_temp/06', 'image');
DTI(:,:,6) = image;
load('~hahntobi/matlab/work_images/DTI_temp/07', 'image');
DTI(:,:,7) = image;

load('~hahntobi/matlab/work_images/DCE_temp/01', 'image');
DCE(:,:,1) = image;
load('~hahntobi/matlab/work_images/DCE_temp/02', 'image');
DCE(:,:,2) = image;
load('~hahntobi/matlab/work_images/DCE_temp/03', 'image');
DCE(:,:,3) = image;
load('~hahntobi/matlab/work_images/DCE_temp/04', 'image');
DCE(:,:,4) = image;
load('~hahntobi/matlab/work_images/DCE_temp/05', 'image');
DCE(:,:,5) = image;
load('~hahntobi/matlab/work_images/DCE_temp/06', 'image');
DCE(:,:,6) = image;

DCE_image = DCE(:,:,1);

DTI_avg = 0;
for i = 2:7
    DTI_avg = DTI(:,:,i)+DTI_avg;
end
DTI_avg = DTI_avg/6;

DTI_avg_large = zeros(512,512);
% enlarge DTI image
for i = 1:256
    for j = 1: 256
        DTI_avg_large(2*i,2*j) = DTI_avg(i,j);
        DTI_avg_large(2*i-1,2*j) = DTI_avg(i,j);
        DTI_avg_large(2*i-1,2*j-1) = DTI_avg(i,j);
        DTI_avg_large(2*i,2*j-1) = DTI_avg(i,j);
    end
end
%figure; imshow(DTI_avg_large,[]); title('DTI average image enlarged');

DTI_avg = DTI_avg_large;
save('/home/hahntobi/matlab/work_images/dti_avg', 'DTI_avg_large');
save_filename = ['/home/hahntobi/matlab/work_images/dti_avg', '.tiff'];
imwrite(DTI_avg_large,save_filename,'tif')

DTI_edge = edge(DTI_avg); % find the edges of the image
%figure; imshow(DTI_edge); title('DTI image edges');

% Insert more pixels in the ADC image to make sure that there has been no
% information generated "out of nothing", since the resolution of the DTI
% image originally has only been 256*256!
DTI_edge2 = DTI_edge;
for i = 1:512
    for j=1:512
        if DTI_edge2(i,j)==1
            DTI_edge2(i+mod(i,2),j+mod(j,2)-1)=1;
            DTI_edge2(i+mod(i,2),j+mod(j,2))=1;
            DTI_edge2(i+mod(i,2)-1,j+mod(j,2))=1;
            DTI_edge2(i+mod(i,2)-1,j+mod(j,2)-1)=1;
        end
    end
end


% Shifting the DTI image in vertical direction to see if there is a better
% fit with the DCE image, which would be an indication for an overall shift
% in vertical direction.
% DTI_edge_3D = zeros(size(DTI_edge,1),size(DTI_edge,2),4);
% DTI_edge_3D(:,:,1)=DTI_edge2;
% for k = 2: 5
%     DTI_edge_3D(:,:,k) = zeros(size(DTI_edge,1),size(DTI_edge,2));
%     for i = 1:512
%         for j=1:512
%             if i~=1
%                 DTI_edge_3D(i,j,k)=DTI_edge_3D(i-1,j,k-1);
%             else
%                 DTI_edge_3D(i,j,k)=DTI_edge_3D(512,j,k-1);
%             end
%         end
%     end 
% end


DCE_edge = edge(DCE_image);
%figure; imshow(DCE_edge,[]); title('DCE image edges');

%figure; imshow(DCE_edge+2*DTI_edge,[]); colormap(jet)
figure; imshow(DCE_edge+2*DTI_edge2,[]); colormap(jet) % Figure 1
% for i=1:4
%     figure; imshow(DCE_edge+2*DTI_edge_3D(:,:,i+1),[]); colormap(jet)
% end

figure; imshow(DCE_edge + 2*shift_image(shift_image(DTI_edge2,0,+2),2,1),[]); colormap(jet) % Figure 2
figure; imshow(DCE_edge + 2*shift_image(shift_image(DTI_edge2,2,-2),2,1),[]); colormap(jet) % Figure 3
figure; imshow(DCE_edge + 2*shift_image(shift_image(DTI_edge2,2,+2),2,1),[]); colormap(jet) % Figure 4
figure; imshow(DCE_edge + 2*shift_image(shift_image(DTI_edge2,1,+2),2,1),[]); colormap(jet) % Figure 5
