function afni_preprocess_afni(flist,anat,f_volreg)
% afni_preprocess_afni(flist[,anat,f_volreg])
% data preprocessing in afni.
% flist: folders of the functional data, in the same sequence as in
% acquisition, excluding '+orig'.
% nslices: number of slices of the functional scans. default: 20.
% TR: unit ms. Default 2000.
% f_volreg: the ind for the file in flist to be used as base for motion
% correction.  The last image in that file will be used as the base.
% default: the last one.
% anat: folder for anatomical data, excluding '+orig'.

tic;

flist = str2cell(flist);
if ~exist('f_volreg','var')
    f_volreg = length(flist);
end

no_anat = false;
if ~exist('anat','var') || isempty(anat)
    no_anat = true;
end
    
nf = length(flist);
[ERR,info] = BrikInfo([flist{nf},'+orig']);
nTR_last = info.DATASET_RANK(2);

for i=1:nf
    cmd = sprintf('3dTshift -rlt+ -prefix s%d_tsh %s+orig',i,flist{i});
    unix(cmd);
end    
    
% volume register
for i=1:nf
   cmd = sprintf('3dvolreg -base s%d_tsh+orig''[%d]'' -prefix s%d_vr -1Dfile motion_s%d.1D s%d_tsh+orig',f_volreg,nTR_last-1,i,i,i); 
   unix(cmd);
end

% detrend 
for i=1:nf
   cmd = sprintf('3dDetrend -prefix s%d_vr_dtr -polort 2 s%d_vr+orig',i,i);
   unix(cmd);
end

%generate brainmask.
cmd = sprintf('3dAutomask -prefix brainmask s%d_vr+orig',nf);
unix(cmd);

% normalized signal intensity
for i=1:nf
    cmd = sprintf('3dClipLevel -mfrac 0.55 s%d_vr+orig',i);
    [s,clip] = unix(cmd);
    clip=str2double(clip);
    disp(clip);
    cmd =  sprintf('3dTstat -prefix s%d_mean s%d_vr+orig',i,i);
    unix(cmd); 
    cmd = sprintf('3dcalc -a s%d_vr_dtr+orig -b s%d_mean+orig -expr ''a/b*step(b-%f)'' -prefix s%d_dtr_norm',i,i,clip,i);
    unix(cmd);
   
end

% align anatomy to functional
if ~no_anat
cmd = sprintf('align_epi_anat.py -epi s%d_vr+orig -anat %s+orig -epi_base %d -partial_axial',f_volreg,anat,nTR_last-1);
unix(cmd);
end

fprintf('Preprocessing finish in %3.1f minutes\n',toc/60);

function re_order(folder)
    disp(['re-order ',folder]);
    old_dir =  cd(folder);
    dir_str = dir('*MRDC*');
    
    dir_str_test = dir('i*MRDC*');
    if isempty(dir_str_test)
        disp('files already ordered');
        return;
    end
    cmd = sprintf('re-order %d',length(dir_str));
    unix(cmd);
    cd(old_dir);
    
    
    
    