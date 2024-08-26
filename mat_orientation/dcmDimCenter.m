function [vsize,center,orient,pos,rotmat]=dcmDimCenter(pos,OrientPat,voxsize,columns,rows)
%dcmDimCenter(dname) or dcmDimCenter(pos,ImageOrientPatient,voxsize)
% pos: 3*nslice; 
% orient: 1*6; ImageOrientationPatient; columns then rows, then slice
% voxsize: 1*2; pixsize from dicom header; for dimension along rows (2nd dim)
% and then along columns (1st dim)
%

% output:
% vsize: 1*3; rows then columns; same for orient label 
% center: 1*3
% orient: label
% pos: 3*1; position of pixel 1,1,1
% rotmat: 3*3 matrix; physical position is rotmat*([i-1,j-1,k-1].*voxsize)'+pos; dicom convention;
% % center of the image: for the first two dimensions the center is on the
% n/2+1 (even n) or (n+1)/2 (odd n) pixels, and for third dimension on
% (n+1)/2, regardless of n;  This is the same as the Siemens scanner
% convention and true for both 2D and 3D sequences.
% dicom convention: 
%  x=R-L, y=A-P, and z=I-S (R,A,I are < 0; L,P,S are > 0).
% nifti convention: L-R; P-A; I-S; 



% center=get_center(pos,voxsize,OrientPat,columns,rows);
% 
% orient=orientLabel(reshape(vori,[3,2])',pos);
%  
warning off;

if isa(pos,'char')
    [pos,OrientPat,voxsize,columns,rows]=get_pos(pos);
end
        
%%
center=get_center(pos,voxsize,OrientPat,columns,rows);

vsize=voxsize;
if size(pos,2)>1
vsize(3)=sos(pos(:,1)-pos(:,2),1);
else 
    vsize(3)=1;%arbitrary
end
rotmat=reshape(OrientPat,[3,2]);
if size(pos,2)>1
    rotmat(:,3)=(pos(:,2)-pos(:,1))/vsize(3);
else
    rotmat(:,3)=cross(rotmat(:,2),rotmat(:,1)); % arbitrary % swithced 1 and 2 on 7/24/2020 to be consistent with slice direction in Siemens .pro file
end

orient=orientLabel(rotmat);
%% switch 1 and 2 [vsize,center,orient,pos,rotmat]; because in dicom column is the first dim but in matlab row is the first dim
vsize=vsize([2,1,3]);
%center=center([2,1,3]);
orient=orient([2,1,3]);
%pos=pos(:,[2,1,3]);
rotmat=rotmat(:,[2,1,3]);
pos=pos(:,1);




function ori=orientLabel(v)

label='LRPASI';  % the negative side


for i=1:size(v,2)
   
    [tmp,ind]=max(abs(v(:,i)));

    j=ind*2+(sign(v(ind,i))-1)/2;
    ori(i)=label(j);
%     f2(i)=ind;
    
end

% 
% a=setdiff(1:3,f2);
% if length(pos)>1
%     dpos=pos(2)-pos(1);
% else
%     dpos=1;  %does not matter
% end
% 
% j=a*2+(sign(dpos)-1)/2;
% ori(3)=label(j);
% 


function [pos,vori,vsize,columns,rows]=get_pos(dname)

% pos: 3*n
% vsize: 1*2
% vroi: 1*6
% The first element for colum
dir_str=dir2([dname,'/*']);
   
pos=[];
nfiles2read=length(dir_str);

i=1;
i2=1;
while i2<=nfiles2read
   if mod(i,10)==0
     %  disp(i);
   end
    
   if ~exist('ns','var')
       try
           h=dicominfo(fullfile(dname,dir_str(i2).name));
         
       catch
           i2=i2+1;
           continue;
       end
       
        tmp=h.ImagePositionPatient; % 1*3; the third component is slice location; the first two are location of the top left voxel.
        vori=h.ImageOrientationPatient;
       
        if isempty(pos)
            pos(:,i)=tmp(:); 
        elseif any(tmp(1)~=pos(1,:))  || any(tmp(2)~=pos(2,:)) || any(tmp(3)~=pos(3,:))
            pos(:,i)=tmp(:);  
        else

            ns=length(pos);  % repeated slice locations found; break;        
            break;
        end

  
        if i==1
            dinfo = dicominfo(fullfile(dname,dir_str(i2).name));
            
            vsize(1:2)=dinfo.PixelSpacing;
            vsize(3)=dinfo.SliceThickness;
           % vsize(4)=1;
            rows=dinfo.Rows;% first dimension when read with ri.m
            columns=dinfo.Columns; %second dimension when read with ri.m
            
        end
        
    %    vori=reshape(vori,[3,2]);
       
        
   end
    i=i+1;    
    i2=i2+1;
end


function center=get_center(pos,vsize,vori,columns,rows)

vori=reshape(vori,[3,2]);
vsize=vsize(:);
for i=1:size(pos,2)
    ij_c=floor(double([columns,rows])/2);        
    center(:,i)=pos(:,i)+(vori.*repmat(vsize(1:2)',[3,1]))*ij_c';  
end

center=mean(center,2);

  
  
  



