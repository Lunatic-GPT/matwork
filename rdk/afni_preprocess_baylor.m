function afni_preprocess_baylor(flist,anat,nslices,TR,f_volreg,TR_volreg)
% afni_preprocess(flist[,anat,nslices,TR,f_volreg,TR_volreg])
% data preprocessing in afni.
% flist: folders of the functional data, in the same sequence as in acquisition  
% nslices: number of slices of the functional scans. default: 20.
% TR: unit ms. Default 2000.
% f_volreg: the ind for the file in flist to be used as base for motion
% correction.  The last image in that file will be used as the base.
% default: the last one.
% anat: folder for anatomical data.

tic;
if ~exist('nslices','var')
    nslices=29;
end

if ~exist('TR','var')
    TR = 2000;
end

if ~exist('f_volreg','var')
    f_volreg = 1;
end

if ~exist('TR_volreg','var')
    TR_volreg = 0;
end

no_anat = false;
if ~exist('anat','var') || isempty(anat)
    no_anat = true;
end
    
afni_dir = pwd;

% import data and slice timing correction
cd(afni_dir);

if ~no_anat
  cmd = sprintf('to3d -anat -prefix anat %s',fullfile('..',anat,'*'));
  unix(cmd);
end

flist = str2cell(flist);
nf = length(flist);
for i=1:nf
    dir_str = dir(fullfile('..',flist{i},'*.ima'));
    nTR = length(dir_str);
    
    if mod(nslices,2) == 1
        cmd = sprintf('to3d -epan -assume_dicom_mosaic -time:zt %d %d %d alt+z -prefix s%d %s',nslices,nTR,TR,i,fullfile('..',flist{i},'*.ima'));
    else
        cmd = sprintf('to3d -epan -assume_dicom_mosaic -time:zt %d %d %d alt+z2 -prefix s%d %s',nslices,nTR,TR,i,fullfile('..',flist{i},'*.ima'));
    end
    if ~exist(sprintf('s%d+orig.HEAD',i),'file')
     unix(cmd);
    end
    cmd = sprintf('3dTshift -rlt+ -tzero 0 -prefix s%d_tsh s%d+orig',i,i);  % aligned to the beginning of each TR 04/21/2010, X.Z.
    if ~exist(sprintf('s%d_tsh+orig.HEAD',i),'file')
     unix(cmd);
    end
end    
    

% volume register
for i=1:nf
   cmd = sprintf('3dvolreg -base s%d_tsh+orig''[%d]'' -prefix s%d_vr -1Dfile motion_s%d.1D s%d_tsh+orig',f_volreg,TR_volreg,i,i,i); 
   if ~exist(sprintf('s%d_vr+orig.HEAD',i),'file')
     unix(cmd);
    end
   cmd = sprintf('3dmerge -1blur_fwhm 3.4 -doall -prefix s%d_blur s%d_vr+orig',i,i);
   if ~exist(sprintf('s%d_blur+orig.HEAD',i),'file')
     unix(cmd);
    end
end

% detrend 
for i=1:nf
  % cmd = sprintf('3dDetrend -prefix s%d_vr_dtr -polort 2 s%d_vr+orig',i,i);
   cmd = sprintf('3dDetrend -prefix s%d_vr_dtr -polort 2 s%d_blur+orig',i,i);
   if ~exist(sprintf('s%d_vr_dtr+orig.HEAD',i),'file')
     unix(cmd);
   end
end


% normalized signal intensity
clipv = zeros(1,nf);
for i=1:nf
    cmd = sprintf('3dClipLevel -mfrac 0.55 s%d_blur+orig',i);
    [s,clip] = unix(cmd);
    clipv(i)=str2double(clip(end-4:end));
    disp(clip);
    cmd =  sprintf('3dTstat -prefix s%d_mean s%d_blur+orig',i,i);
    if ~exist(sprintf('s%d_mean+orig.HEAD',i),'file')
     unix(cmd);
    end 
    cmd = sprintf('3dcalc -a s%d_vr_dtr+orig -b s%d_mean+orig -expr ''a/b*step(b-%f)'' -prefix s%d_dtr_norm',i,i,clipv(i),i);
    if ~exist(sprintf('s%d_dtr_norm+orig.HEAD',i),'file')
     unix(cmd);
    end
   
end

%generate brainmask.
dsets = '';
expr = sprintf('step(a-%f)',clipv(1));
for i=1:nf
  dsets = sprintf('%s -%c s%i_mean+orig',dsets,'a'+i-1,i);
  if i>1
    expr = sprintf('%s*step(%c-%f)',expr,'a'+i-1,clipv(i));
  end
end
cmd = sprintf('3dcalc %s -expr ''%s'' -prefix brainmask',dsets,expr);
if ~exist('brainmask+orig.HEAD','file')
     unix(cmd);
end

% align anatomy to functional
if ~no_anat
cmd = sprintf('align_epi_anat.py -epi s%d_vr+orig -anat anat+orig -epi_base %d -partial_axial',f_volreg,TR_volreg);
%unix(cmd);
end

fprintf('Preprocessing finish in %3.1f minutes\n',toc/60);

    
    
    
    
