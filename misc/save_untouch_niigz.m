function save_untouch_niigz(nii,prefix,newdata)
%save_untouch_niigz(nii,prefix,savemat)
% 
% if ~exist('savemat','var')
%     savemat=false;
% end

if isa(nii,'char') %file name
    nii=load_untouch_niigz(nii);
end

if exist('newdata','var')
    nii.img=newdata;
end

if length(prefix)>4 && strcmp(prefix(end-3:end),'.nii')
  prefix=prefix(1:end-4);
end

if length(prefix)>7 && strcmp(prefix(end-6:end),'.nii.gz')
  prefix=prefix(1:end-7);
end

if isa(nii.img,'int16')
    nii.hdr.dime.datatype=4;
    nii.hdr.dime.bitpix=16;
elseif isa(nii.img,'single')
    nii.hdr.dime.datatype=16;
    nii.hdr.dime.bitpix=32;
elseif isa(nii.img,'double')
    nii.hdr.dime.datatype=64;
    nii.hdr.dime.bitpix=64;
elseif isa(nii.img,'bool')
    nii.hdr.dime.datatype=1;
    nii.hdr.dime.bitpix=1;
end
 
nii.hdr.dime.dim(5)=size(nii.img,4);
nii.hdr.dime.dim(1)=ndims(nii.img);
nii.hdr.dime.dim(2)=size(nii.img,1);
nii.hdr.dime.dim(3)=size(nii.img,2);
nii.hdr.dime.dim(4)=size(nii.img,3);

save_untouch_nii(nii,[prefix,'.nii']);
gzip([prefix,'.nii']);
delete([prefix,'.nii']);
% 
% if savemat
%     d=nii;
%   save([prefix,'.mat'],'d');
% end


