function afni_preprocess_pc_rvtc(flist,anat,resp,ecg,nslices,TR,f_volreg)
% afni_preprocess_pc_rvtc(flist[,anat,resp,ecg,nslices,TR,f_volreg])
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
    nslices=20;
end

if ~exist('TR','var')
    TR = 2000;
end

if ~exist('f_volreg','var')
    f_volreg = length(flist);
end

no_anat = false;
if ~exist('anat','var') || isempty(anat)
    no_anat = true;
end
    
afni_dir = pwd;

if ~no_anat
  disp('re-ordering anatomical files');
  re_order(fullfile('..',anat));
end

flist=str2cell(flist);
nf = length(flist);
disp('re-ordering functional files');
for i= 1:nf
   re_order(fullfile('..',flist{i}));
end

% import data and slice timing correction
cd(afni_dir);

if ~no_anat
  cmd = sprintf('to3d -anat -prefix anat %s',fullfile('..',anat,'*'));
   if ~exist('anat+orig.HEAD','file')
    unix(cmd);
   end
  
end

for i=1:nf
    
    if ~exist(sprintf('s%d+orig.HEAD',i),'file')
       cmd = sprintf('to3d -epan -time:zt %d %d %d alt+z -prefix s%d %s',nslices,ns/nslices,TR,i,fullfile('..',flist{i},'MRDC*'));
       unix(cmd);
    end
    
    if ~exist(sprintf('pcs%d+orig.HEAD',i),'file')
     retroicor(sprintf('s%d+orig.HEAD',i),resp{i},ecg{i},2)
    end
    
    if ~exist(sprintf('rvtc_pcs%d+orig.HEAD',i),'file')
     rvtc(sprintf('pcs%d+orig.HEAD',i),resp{i})
    end
  %  cmd = sprintf('3dTshift -rlt+ -prefix s%d_tsh s%d+orig',i,i);
    %cmd = sprintf('3dcopy s%d+orig s%d_tsh',i,i); % no correction
    if ~exist(sprintf('rvtc_pcs%d_tsh+orig.HEAD',i),'file')
      cmd = sprintf('3dTshift -rlt+ -tzero 0 -prefix rvtc_pcs%d_tsh rvtc_pcs%d+orig',i,i);  % aligned to the beginning of each TR 04/21/2010, X.Z.
      unix(cmd);
    end
    
end    
    
[tmp,info]= BrikLoad(sprintf('s%d+orig',nf));
   nTR_last = info.DATASET_RANK(2);

% volume register
for i=1:nf
   if ~exist(sprintf('rvtc_pcs%d_vr+orig.HEAD',i),'file')
    cmd = sprintf('3dvolreg -base rvtc_pcs%d_tsh+orig''[%d]'' -prefix rvtc_pcs%d_vr -1Dfile rvtc_pcmotion_s%d.1D rvtc_pcs%d_tsh+orig',f_volreg,nTR_last-1,i,i,i);
    unix(cmd);
   end
   
end

% detrend 
for i=1:nf
    if ~exist(sprintf('rvtc_pcs%d_vr_dtr+orig.HEAD',i),'file')
     cmd = sprintf('3dDetrend -prefix rvtc_pcs%d_vr_dtr -polort 2 rvtc_pcs%d_vr+orig',i,i); 
     unix(cmd);
   end
end


% normalized signal intensity
clipv = zeros(1,nf);
for i=1:nf
    
    if ~exist(sprintf('rvtc_pcs%d_mean+orig.HEAD',i),'file')
     cmd =  sprintf('3dTstat -prefix rvtc_pcs%d_mean rvtc_pcs%d_vr+orig',i,i);
     unix(cmd); 
    end
    
    cmd = sprintf('3dClipLevel -mfrac 0.55 rvtc_pcs%d_vr+orig',i);
        [s,clip] = unix(cmd);
        clipv(i)=str2double(clip);
        disp(clip);
        
    if ~exist(sprintf('rvtc_pcs%d_dtr_norm+orig.HEAD',i),'file') 
        cmd = sprintf('3dcalc -a rvtc_pcs%d_vr_dtr+orig -b rvtc_pcs%d_mean+orig -expr ''a/b*step(b-%f)'' -prefix rvtc_pcs%d_dtr_norm',i,i,clipv(i),i);
        unix(cmd);
    end
   
end

%generate brainmask.
if ~exist('brainmask+orig.HEAD','file')
dsets = '';
expr = sprintf('step(a-%f)',clipv(1));
for i=1:nf
  dsets = sprintf('%s -%c rvtc_pcs%i_mean+orig',dsets,'a'+i-1,i);
  if i>1
    expr = sprintf('%s*step(%c-%f)',expr,'a'+i-1,clipv(i));
  end
end
cmd = sprintf('3dcalc %s -expr ''%s'' -prefix brainmask',dsets,expr);
unix(cmd);
end

% align anatomy to functional
if ~no_anat && ~exist('anat_al+orig.HEAD','file')
cmd = sprintf('align_epi_anat.py -epi rvtc_pcs%d_vr+orig -anat anat+orig -epi_base %d -partial_axial',f_volreg,nTR_last-1);
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
    
    
    
    
