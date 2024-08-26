function res=read_boxroi_amide(fname,uname,vsz_new,dim_new)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

xml2=sprintf('roi_%s.xml',fname);
os=readp(xml2,'coordinate_space_offset');
boxsz=readp(xml2,'corner');

xmlu=['data-set_',uname,'_raw-data.xml'];
xmlu2=['data-set_',uname,'.xml'];

dim2=readp(xmlu,'dim');
os2=readp(xmlu2,'coordinate_space_offset');
vsz=readp(xmlu2,'voxel_size');


if exist('dim_new','var')
    dim2=dim_new;
    vszu=vsz_new;
end

sz2=readp(xmlu2,'coordinate_space_z');
sy2=readp(xmlu2,'coordinate_space_y');
sx2=readp(xmlu2,'coordinate_space_x');

if sx2(1)==-1
    os2(1)=-os2(1);
end

if sy2(2)==-1
    os2(2)=-os2(2);
end

if sz2(3)==-1
    os2(3)=-os2(3);
end

dim=ceil(boxsz./vsz);
xi=round((os-os2)./vsz+1);


res=zeros(dim2(1:3));
res(xi(1):xi(1)+dim(1)-1,xi(2):xi(2)+dim(2)-1,xi(3):xi(3)+dim(3)-1)=1;

roi=res;
if ~exist([fname,'.mat'],'file')
  save(fname,'roi');
end

%res=flipdim(res,2);

% 
% for i=1:size(res,3)
%     if any(vec(res(:,:,i))>0)
%         roi=res(:,:,i);
%         save(sprintf('%sX%d.mat',fname,i),'roi');
%     end
% end


       

function res=readp(fname,par)

d=xmlread(fname);
thisListItem = d.getElementsByTagName(par);
childNode = thisListItem.item(0).getFirstChild;

res=str2num(childNode.getNodeValue);

if isempty(res)
    res=childNode.getNodeValue;
end

