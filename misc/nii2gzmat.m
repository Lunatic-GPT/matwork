function nii2gzmat(savemat)

if ~exist('savemat','var')
    savemat=1;
end

dir_str=dir('*.nii');

for i=1:length(dir_str)
    
    prefix=strtok(dir_str(i).name,'.');
nii=load_untouch_nii(dir_str(i).name);

if ~exist([dir_str(i).name,'.gz'],'file')
  gzip(dir_str(i).name);
end

delete(dir_str(i).name);

if savemat
  d=nii;
  save([prefix,'.mat'],'d');
end
end