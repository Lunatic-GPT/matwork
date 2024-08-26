function show_wash_in_out_callback(params)


nslices = get(params,'Total number of slices');
fslice =get(params,'First slice');
folder = get(params,'Select image directory');  % folder for downloaded images. Also the folder for saving the shifted images.
data_dir = get(params,'Select data directory'); 



for i=fslice:fslice+nslices-1
    
    
h=figure;
set(h,'Units','pixels');
set(h,'Position',[114,727,1314,298]);

fname = sprintf('%02d_2.mat',i);
temp=load(fullfile(folder,fname),'image');


dir_str = dir(fullfile(data_dir,['*X',num2str(i)]));
temp2 = load(fullfile(data_dir,dir_str.name,'ROIs.mat'),'rect_pos','roi','tissue_BW_large','tissue2_BW_large','tissue2_BW_large_b');

rect_pos = temp2.tissue2_BW_large_b;
c_image = full2crop(temp.image,rect_pos);
roi = crop2full(temp2.roi,temp2.rect_pos,temp2.tissue_BW_large);
roi = full2crop(roi,rect_pos);

%tissue = full2crop(temp2.tissue_BW_large,rect_pos);
%tissue2 = full2crop(temp2.tissue2_BW_large,rect_pos);


fname = sprintf('%02d_wash_in_norm_diff.mat',i);
wi = load(fullfile(folder, fname),'image');
wi = full2crop(wi.image,rect_pos);

fname = sprintf('%02d_wash_out.mat',i);
wo = load(fullfile(folder, fname),'image');
wo = full2crop(wo.image,rect_pos);


in = round(max(c_image(:)));
subplot(1,4,1),subimage(c_image,gray(in));
 title(['Slice ',num2str(i),', t=2']);
axis off image;
 set(gca,'Position',[0,0,0.2,0.9]);

norm_image = 100*c_image/max(c_image(:));
norm_image(norm_image<2) = 2;
%edge1 = bwmorph(roi,'remove') | bwmorph(imfill(tissue,'hole'),'remove') | bwmorph(imfill(tissue2,'hole'),'remove');
edge1 = bwmorph(roi,'remove');

norm_image(edge1>0) =1;
cm = gray(100);
cm(1,:) = [1,0,0];
subplot(1,4,2),subimage(norm_image,cm);
title(['Slice ',num2str(i),' with roi']);
axis off;
set(gca,'Position',[0.25,0,0.2,0.9]);

    jet_inv = flipud(jet(90));   
    colormap(jet_inv);
    wi_min = min(wi(:));
    wi_max=max(wi(:));
    
    subplot(1,4,3), subimage((wi-wi_min)*90/(wi_max-wi_min),jet_inv);
    title(['Slice ',num2str(i),': Wash-in (normalized intensity increase)']);
    tlabel = {};
    for it=1:4
        tlabel{it} =sprintf('%3.2f',(it-1)*(wi_max-wi_min)/3+wi_min);
    end
    colorbar('YTickLabel',tlabel,'YTick',1:89/3:90);
    axis off;
    set(gca,'Position',[0.5,0,0.20,0.9]);
    
    jet_inv = flipud(jet(180));
    subplot(1,4,4), subimage(wo+91,jet_inv);
    title(['Slice ',num2str(i),': Wash-out']);
    colorbar('YTickLabel',{'-90','-45','0','45','90'},'YTick',1:89/4:90);
    set(gca,'Position',[0.75,0,0.2,0.9]);
    axis off;
    
end


