function [tissue_BW,tissue2_BW]=tissue(c_image,roi,folder)
%% Create Tissue Area if we are satisfied with the ROI


BW_morph_inner = bwmorph(roi,'remove'); % we are only interested in the margin pixels


%% TEST-GAP
lesion_area = sum(roi(:)>0);
fil = ones(3,3);
TG1 = filter2(fil,BW_morph_inner);
TG1(TG1>0)=1;


%% Only a try!!!!

    gap_BW = imfill(TG1, 'holes');
    gap_copy = gap_BW-roi; % THE GAP

tissue_morph = gap_BW;
tissue_area = 0;
size_flag = 0;
ratio = 0;
while size_flag == 0
    TG1 = filter2(fil,tissue_morph);
    TG1(TG1>0)=1;
    tissue_morph_old = tissue_morph;
    %tissue_morph = bwmorph(TG1,'remove');
    tissue_morph = TG1;
    tissue_area_old = tissue_area;
    tissue_copy = tissue_morph - gap_BW; 
    tissue_area = sum(sum(tissue_copy==1));
    
    if tissue_area > lesion_area
        outer_layer = bwmorph(tissue_morph,'remove');
        BW_label_morph = bwlabel(outer_layer);
        pixels = regionprops(BW_label_morph, 'PixelList');
        if (tissue_area - lesion_area) > (lesion_area - tissue_area_old) % the old tissue layer is closer to the lesion area than the larger tissue layer
            tissue_copy = tissue_morph_old - gap_BW;
            ratio = lesion_area/sum(sum(tissue_copy==1));
            while ratio ~= 1
                listsize = size(pixels.PixelList,1);
                random = rand;
                random = ceil(random*listsize);
                tissue_morph_old(pixels.PixelList(random,2), pixels.PixelList(random,1))=1; % add pixels
                tissue_copy = tissue_morph_old - gap_BW;
                ratio = lesion_area/sum(sum(tissue_copy==1));
                tissue_morph = tissue_morph_old;
            end
        else
            while ratio ~= 1
                random = ceil(rand*size(pixels.PixelList,1));
                tissue_morph(pixels.PixelList(random,2), pixels.PixelList(random,1))=0; % take off pixels
                tissue_copy = tissue_morph - gap_BW;
                ratio = lesion_area/sum(sum(tissue_copy==1));
            end
        end
        size_flag = 1;
    end
    
    if tissue_area == tissue_area_old
         errordlg('Please use a larger cropped image!  Tissue around lesion is bigger than your chosen cropped image.');
         error('Please use a larger cropped image!  Tissue around lesion is bigger than your chosen cropped image.');
    end
end



%%

tissue_BW = imfill(tissue_morph,'holes');
tissue_BW_filled = tissue_BW; % copy for tissue 2 calculation
tissue_BW = tissue_BW - gap_BW;

%% NEW (4/1/08): TISSUE 2
%tissue_morph = bwmorph(gap_BW,'remove'); % we are only interested in the margin pixels
tissue2_morph = tissue_BW_filled;
tissue2_area = 0;
size_flag = 0;
ratio = 0;
while size_flag == 0
    TG1 = filter2(fil,tissue2_morph);
    TG1(TG1>0)=1;
    tissue2_morph_old = tissue2_morph;
    %tissue2_morph = bwmorph(TG1,'remove');
    tissue2_morph = TG1;
    tissue2_area_old = tissue2_area;
    tissue2_copy = tissue2_morph - tissue_BW_filled; % only the tissue area 
    tissue2_area = sum(sum(tissue2_copy==1));
    
    if tissue2_area > lesion_area
        outer_layer = bwmorph(tissue2_morph,'remove');
        BW_label_morph = bwlabel(outer_layer);
        pixels = regionprops(BW_label_morph, 'PixelList');
        if (tissue2_area - lesion_area) > (lesion_area - tissue2_area_old) % the old tissue2 layer is closer to the lesion area than the larger tissue2 layer
            tissue2_copy = tissue2_morph_old - tissue_BW_filled;
            ratio = lesion_area/sum(sum(tissue2_copy==1));
            while ratio ~= 1
                listsize = size(pixels.PixelList,1);
                random = rand;
                random = ceil(random*listsize);
                tissue2_morph_old(pixels.PixelList(random,2), pixels.PixelList(random,1))=1; % add pixels
                tissue2_copy = tissue2_morph_old - tissue_BW_filled;
                ratio = lesion_area/sum(sum(tissue2_copy==1));
                tissue2_morph = tissue2_morph_old;
            end
        else
            while ratio ~= 1
                random = ceil(rand*size(pixels.PixelList,1));
                tissue2_morph(pixels.PixelList(random,2), pixels.PixelList(random,1))=0; % take off pixels
                tissue2_copy = tissue2_morph - tissue_BW_filled;
                ratio = lesion_area/sum(sum(tissue2_copy==1));
            end
        end
        size_flag = 1;
    end
    
    if tissue2_area == tissue2_area_old
         errordlg('Please use a larger cropped image!  Tissue2 around lesion is bigger than your chosen cropped image.');
         error('Please use a larger cropped image!  Tissue2 around lesion is bigger than your chosen cropped image.');
    end
end


tissue2_BW = imfill(tissue2_morph,'holes');
tissue2_BW = tissue2_BW - tissue_BW_filled;
%% END TISSUE 2


STATS = regionprops(tissue2_BW, 'BoundingBox');
temp = STATS.BoundingBox;

x_1 = floor(temp(2));
x_2 =x_1+ ceil(temp(4));
y_1 = floor(temp(1));
y_2 = ceil(temp(3))+y_1;

if (x_1 == 0) || (y_1 == 0) || (x_2 >size(c_image,1)) || (y_2>size(c_image,2))
    errordlg('Please select a larger rectangle around the lesion!  Tissue2 around lesion is bigger than your chosen cropped image.');
end


if exist('folder','var')

inner_BW = roi;
filename = fullfile(folder, 'tissue_BWs');
save(filename, 'tissue_BW', 'inner_BW', 'tissue_copy', 'tissue_morph', 'gap_BW','gap_copy', 'tissue2_BW', 'tissue2_copy', 'tissue2_morph');


end
