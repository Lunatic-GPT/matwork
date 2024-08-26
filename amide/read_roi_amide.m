function res=read_roi_amide(fname,uname,vsz_new,dim_new)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

xml1=['roi_',fname,'_map_data.xml'];
format=readp(xml1,'raw_format');
dim=readp(xml1,'dim');

fid=fopen(sprintf('roi_%s_map_data.dat',fname));

  if strcmp(format, 'ubyte-8-ne')
        d=fread(fid,'uint8');
  
  end

fclose(fid);

d=reshape(d,dim);
xml2=sprintf('roi_%s.xml',fname);
vsz=readp(xml2,'voxel_size');
os=readp(xml2,'coordinate_space_offset');


xmlu=['data-set_',uname,'_raw-data.xml'];
xmlu2=['data-set_',uname,'.xml'];
vszu=readp(xmlu2,'voxel_size');

dim2=readp(xmlu,'dim');
os2=readp(xmlu2,'coordinate_space_offset');


if exist('dim_new','var')
    dim2=dim_new;
    vszu=vsz_new;
end
res=zeros(dim2(1:3));


%xi=round((os-os2)./vszu+1);

%res(xi(1):xi(1)+dim(1)-1,xi(2):xi(2)+dim(2)-1,xi(3):xi(3)+dim(3)-1)=d;

for i=1:size(d,1)
    for j=1:size(d,2)
        for k=1:size(d,3)
            
            xi=floor((os+[i-1,j-1,k-1].*vsz-os2)./vszu)+1;
            xf=ceil((os+[i,j,k].*vsz-os2)./vszu)+1;
            
            res(xi(1):xf(1),xi(2):xf(2),xi(3):xf(3))=repmat(d(i,j,k),size(xf)-size(xi)+1) | res(xi(1):xf(1),xi(2):xf(2),xi(3):xf(3));
            
        end
    end
end
%res=flipdim(res,2);

  roi=res;
  if ~exist(fname,'file')
      
   save([fname,'.mat'],'roi');
  end

        
        





function res=readp(fname,par)

d=xmlread(fname);
thisListItem = d.getElementsByTagName(par);
childNode = thisListItem.item(0).getFirstChild;

res=str2num(childNode.getNodeValue);

if isempty(res)
    res=childNode.getNodeValue;
end

