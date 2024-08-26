function to3d_siemens(flist,nslices,TR)
% afni_preprocess(flist[,nslices,TR])
% data preprocessing in afni.
% flist: folders of the functional data, in the same sequence as in acquisition  
% nslices: number of slices of the functional scans. default: 20.
% TR: unit ms. Default 2000.
% f_volreg: the ind for the file in flist to be used as base for motion
% correction.  The last image in that file will be used as the base.
% default: the last one.
% anat: folder for anatomical data.

tic;


flist = str2cell(flist);
nf = length(flist);
for i=1:nf
    
     
    dir_str = dir(fullfile(flist{i},'*.IMA'));
   
    hdr=dicominfo(fullfile(flist{i},dir_str(1).name));
    if ~exist('TR','var')
     TR = hdr.RepetitionTime;
    end

if ~exist('nslices','var')
    nslices=hdr.Private_0019_100a;
end


    nTR = length(dir_str);
    
    if mod(nslices,2) == 1
        cmd = sprintf('to3d -epan -assume_dicom_mosaic -time:zt %d %d %d alt+z -prefix %s %s',nslices,nTR,TR,flist{i},fullfile(flist{i},'*.IMA'));
    else
        cmd = sprintf('to3d -epan -assume_dicom_mosaic -time:zt %d %d %d alt+z2 -prefix %s  %s',nslices,nTR,TR,flist{i},fullfile(flist{i},'*.IMA'));
    end
    
unix(cmd);
end    
    


    
    
    
    
