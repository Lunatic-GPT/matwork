function res = read_data_amide_regrid(fname,uname)
%read data fname and resample it according to the grid defined in uname.

d=read_data_amide(fname);

xml2=sprintf('data-set_%s.xml',fname);
vsz=readp(xml2,'voxel_size');
os=readp(xml2,'coordinate_space_offset');
sz=readp(xml2,'coordinate_space_z');
sy=readp(xml2,'coordinate_space_y');
sx=readp(xml2,'coordinate_space_x');

xmlu=['data-set_',uname,'_raw-data.xml'];
xmlu2=['data-set_',uname,'.xml'];
vszu=readp(xmlu2,'voxel_size');

dim2=readp(xmlu,'dim');
os2=readp(xmlu2,'coordinate_space_offset');
sz2=readp(xmlu2,'coordinate_space_z');
sy2=readp(xmlu2,'coordinate_space_y');
sx2=readp(xmlu2,'coordinate_space_x');

vszu=vszu(1:3).*[sx2(1),sy2(2),sz2(3)];
vsz=vsz.*[sx(1),sy(2),sz(3)];

res=zeros(dim2(1:3));

for i=1:size(res,1)
    for j=1:size(res,2)
        for k=1:size(res,3)
            
            xi=ceil((os2-os+[i-1,j-1,k-1].*vszu)./vsz+0.5);
            xf=floor((os2-os+[i,j,k].*vszu)./vsz+0.5);
            
            tmp=xi;
            xi(xi>xf)=xf(xi>xf);
            xf(xf<xi)=tmp(xf<xi);
            
            xi(xi<1)=1;
            if any(xf<1)
                continue;
            end
            
            %xf(xf>xi)=xi(xf>xi);
            if any(xf>size(d))
                continue;
            end
            
            res(i,j,k)=mean(vec(d(xi(1):xf(1),xi(2):xf(2),xi(3):xf(3))));
            
        end
    end
end
%res=flipdim(res,2);
roi=res;
  fname=sprintf('%s_to_%s.mat',fname,uname);
 save(fname,'roi');






function res=readp(fname,par)

d=xmlread(fname);
thisListItem = d.getElementsByTagName(par);
childNode = thisListItem.item(0).getFirstChild;

res=str2num(childNode.getNodeValue);

if isempty(res)
    res=childNode.getNodeValue;
end

