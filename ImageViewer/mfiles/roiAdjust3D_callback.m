function roiAdjust3D_callback(params)

nscreen = 2;

first_slice = get(params, 'first slice');
nslice = get(params, 'total number of slices');
data_dir = get(params, 'Select data directory name');


dir_struct = dir(fullfile(data_dir,'*X*'));
if isempty(dir_struct)
    error(['Wrong data directory: ', data_dir]);
end
set(gcbo,'Enable','off')
rect_pos_arr = cell(1,nslice);
c_image_arr = cell(1,nslice);
fin_ROI_arr = cell(1,nslice);
for slice = 1:nslice
    dir_struct = dir(fullfile(data_dir,['*X',num2str(slice+first_slice-1)]));
    
    if(length(dir_struct)>1)
        disp(['Multiple directories for slice ',num2str(slice), ' in ',data_dir]);
        disp(['Use ',dir_struct(1).name]);
        dir_struct = dir_struct(1);
    end
    
    [token,rem] = strtok(dir_struct.name, 'X');
    if isempty(rem)
        token=rem;
    end
    DCE_dir = fullfile(data_dir,[token 'X' num2str(slice+first_slice-1)]);
    
    fname = fullfile(DCE_dir,'original_image.mat');
    load(fname,'o_image','c_image');
    
    fname = fullfile(DCE_dir,'ROIs.mat');
    load(fname,'roi','inner','outer','rect_pos');
                    
    if slice ==1
      image3d=zeros([size(o_image),nslice]);
      roi3d = zeros([size(o_image),nslice]);
      inner3d = zeros([size(o_image),nslice]);
      outer3d = zeros([size(o_image),nslice]);
    end
    
    image3d(:,:,slice) = o_image;    
    c_image_arr{slice} =c_image;
    rect_pos_arr{slice} = rect_pos;
    fin_ROI_arr{slice} = roi;
    roi3d(rect_pos(2):(rect_pos(2)+rect_pos(4)),rect_pos(1):(rect_pos(1)+rect_pos(3)),slice)=roi;
    inner3d(rect_pos(2):(rect_pos(2)+rect_pos(4)),rect_pos(1):(rect_pos(1)+rect_pos(3)),slice) = inner;
    outer3d(rect_pos(2):(rect_pos(2)+rect_pos(4)),rect_pos(1):(rect_pos(1)+rect_pos(3)),slice) = outer;
    
end

roi3d = mdetect3D(roi3d,inner3d,outer3d,image3d);
tissue_BW_large = zeros([size(o_image),nslice]);
tissue2_BW_large = zeros([size(o_image),nslice]);

% show both the 2d and 3d rois
for i=1:nslice
   
    modi = mod(i-1,12)+1;
    fleft = (mod(modi-1,3))/3/nscreen;
    flower = 1-ceil(modi/3)/4;
   
  h= figure('Units','Normalized','Position',[fleft,flower,1/3/nscreen,0.20],'MenuBar','none');
    subplot(1,2,1);
       
     image_contour(c_image_arr{i},fin_ROI_arr{i});
     title(['slice ',num2str(i),', 2d roi']);
     axis off image;
     set(gca,'Position',[0,0,0.48,0.9]);
     subplot(1,2,2);
     rect_pos = rect_pos_arr{i};
     roi3d_c = roi3d(rect_pos(2):rect_pos(2)+rect_pos(4),rect_pos(1):(rect_pos(1)+rect_pos(3)),i);
     [tissue_BW,tissue2_BW] =   tissue(c_image_arr{i},roi3d_c);
     
     roi3d_c=adjustROI(roi3d_c,tissue_BW,c_image_arr{i});
     
     tissue_BW_large(:,:,i) = crop2full(tissue_BW,rect_pos,o_image);
     tissue2_BW_large(:,:,i) = crop2full(tissue2_BW,rect_pos,o_image);
     roi3d(:,:,i) = crop2full(roi3d_c,rect_pos,o_image);
     
     
     image_contour(c_image_arr{i},roi3d_c);
     title(['slice ',num2str(i),', 3d roi']);
     axis image off;
     set(gca,'Position',[0.5,0,0.48,0.9]);
     
     if ~exist(fullfile(data_dir,'3dROI'),'dir')
        mkdir( fullfile(data_dir,'3dROI'));
     end
     figname = fullfile(data_dir,'3dROI',['ROI_2d_3d_slice',num2str(i),'.png']);
     saveas(gcf,figname);
     uicontrol('Parent',h,'Style','pushbutton','String','Close','CallBack','close','Units','Normalized','Position',[0.425,0.91,0.15,0.08],'ForegroundColor','red');
     set(h,'Units','Normalized','Position',[fleft,flower,1/3/nscreen,0.23]);
end

 user_entry = questdlg('Do you want to save this 3d ROI?','save ROI?','Yes','No','Yes');

if strcmp(user_entry,'Yes')
    inner_BW_large = roi3d;
    save(fullfile(data_dir,'3dROI','ROI3d.mat'),'inner_BW_large', 'tissue_BW_large', 'tissue2_BW_large','params');
end

set(gcbo,'Enable','on')




