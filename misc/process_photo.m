function process_photo(fname,lr,ud,type)
% fname: image file name
% lr: (1*2) horizontal pixel range of head
% ud: (1*2):vertical pixel range of head
% type: 'passport' or 'td';

% fname='20180526_171301.jpg';
% 
% lr=[457,1049];  %[181,693];  % head range in pixels
% ud=[942,1738];
%  type='td';
 %%
a=imread(fname);
in=imfinfo(fname);
 if strcmp('passport',type) %US passport


pic_dim=[2,2]*25.4;  % mm
photo_dim =[6,4]*25.4; % mm
head_dim_range=[30,35]; %mm vertical dimension range % official [25,35];

 elseif strcmp(type,'travel document')

pic_dim=[48,33];  % mm
photo_dim =[6,4]*25.4; % mm
head_dim_range=[28,33]; %mm vertical dimension range

 else
 pic_dim=[53,35];  % mm
photo_dim =[6,4]*25.4;  % mm
head_dim_range=[28,33]*35/33; %mm vertical dimension range
      
 end

if in.Orientation==6
   a=permute(a,[2,1,3]); 
   a=flip(a,2);
end
gap=10;  % gap between tiled images in pixels
%%
head_dim=mean(head_dim_range,2);

pixel_center=[mean(ud),mean(lr)];
vertical_offset=round(diff(ud)/head_dim(1)*0);  % 0 mm shift of head center from picture center 

npixel=round(diff(ud)/(head_dim(1))*pic_dim);  %number of pixels for the each head pic  

rng1=round(pixel_center-npixel/2)+[vertical_offset,0];
rng2=rng1+npixel;

img=a(rng1(1):rng2(1),rng1(2):rng2(2),:);



img2=255*ones(size(img)+[gap,gap,0]);
img2(gap/2+1:end-gap/2,gap/2+1:end-gap/2,:)=img;

if strcmp('passport',type) 
img3=repmat(img2,[2,1,1]);
elseif strcmp(type,'travel document')
img3=repmat(img2,[2,2,1]);
    
else
   img3=img2;
   
   
end

fname=strtok(fname,'.');
imwrite(uint8(img),[fname,'_',type,'_1.tif'],'TIFF','Resolution',[300,300]);

figure;imshow(img3(:,:,1),[]);
  
photo_pixel=round(photo_dim*npixel(1)/pic_dim(1));
   
sz=size(img3);
img4=uint8(255*ones([photo_pixel,sz(3)]));
  
img4(1:sz(1),1:sz(2),1:sz(3))=uint8(img3);
   
   img5=circshift(img4,400,1);
   img5=circshift(img5,400,2);
   img5(1:10,:,:)=0;
   img5(:,1:10,:)=0;
   img5(end-9:end,:,:)=0;
   img5(:,end-9:end,:)=0;

 
     
   imwrite(img5,[fname,'_',type,'.tif'],'TIFF','Resolution',[300,300]);
   
   
   