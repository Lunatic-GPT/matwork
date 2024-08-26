function roi2=roi2shape(fname,shape)

if isa(fname,'char')
  nii=load_untouch_niigz(fname);
  roi=nii.img;
else
    roi=fname;
end

if strcmp(shape,'circle')
    roi2=roi2circle(roi,3,false);
elseif strcmp(shape,'crosshair')
    roi2=roi2crosshair(roi);
elseif strcmp(shape,'circle_hollow')
    roi2=roi2circle(roi,4,false);
    roi2=roi_boundary(roi2,3);
elseif strcmp(shape,'square')
    roi2=roi2square(roi);
    roi2=roi_boundary(roi2,3);
elseif strcmp(shape,'circle_hollow2')
    roi2=roi2circle(roi,6,false);
    roi2=roi_boundary(roi2,2);
    
end

if isa(fname,'char')
nii.img=roi2;
prefix=[fname(1:end-7),'_',shape];
save_untouch_niigz(nii,prefix);
end

function roi2=roi2crosshair(roi)


roi_orig=roi;
roi=clusterize2(roi>0);

unq=unique(roi(:));
unq(unq==0)=[];
roi2=0*roi;

for i=1:length(unq)
    xy=roiCOM(roi==unq(i));
    m=0*roi;
    m(xy(1)-3:xy(1)+3,xy(2)-1:xy(2)+1)=1;
    m(xy(1)-1:xy(1)+1,xy(2)-3:xy(2)+3)=1; 
        
    roi2(m>0)=mean(roi_orig(roi==unq(i)));
    
end

function roi2=roi2circle(roi,rad,edge_only)

roi_orig=roi;
roi=clusterize2(roi>0);

unq=unique(roi(:));
unq(unq==0)=[];
roi2=0*roi;

for i=1:length(unq)
    center=roiCOM(roi==unq(i));
    m=mask_circle(size(roi),rad,center,1,edge_only);
       
    roi2(m>0)=mean(roi_orig(roi==unq(i)));
end

function roi2=roi2square(roi)

roi_orig=roi;
roi=clusterize2(roi>0);

unq=unique(roi(:));
unq(unq==0)=[];
roi2=0*roi;

for i=1:length(unq)
    xy=roiCOM(roi==unq(i));
    m=0*roi;
    m(xy(1)-4:xy(1)+4,xy(2)-4:xy(2)+4)=1;
    roi2(m>0)=mean(roi_orig(roi==unq(i)));
    
end




