function afni_preprocess_pc(flist,anat,resp,ecg,nslices,TR,name_volreg,i_volreg)
% afni_preprocess_pc(flist[,anat,resp,ecg,nslices,TR,f_volreg])
% data preprocessing in afni.
% flist: folders of the functional data, in the same sequence as in acquisition  
% nslices: number of slices of the functional scans. default: 20.
% TR: unit ms. Default 2000.
% f_volreg: the ind for the file in flist to be used as base for motion
% correction.  The last image in that file will be used as the base.
% default: the last one.
% anat: folder for anatomical data.

flist = str2cell(flist);
resp = str2cell(resp);
ecg = str2cell(ecg);
nf = length(flist);


tic;
if ~exist('nslices','var')
    nslices=20;
end

if ~exist('TR','var')
    TR = 2000;
end

if ~exist('name_volreg','var')
    name_volreg = sprintf('pcs%d_tsh+orig',nf);
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
    dir_str = dir(fullfile('..',flist{i},'*MRDC*'));
    ns = length(dir_str);
    
    if i == nf
        nTR_last = ns/nslices;
    end
    
    cmd = sprintf('to3d -epan -time:zt %d %d %d alt+z -prefix s%d %s',nslices,ns/nslices,TR,i,fullfile('..',flist{i},'MRDC*'));
    
    if ~exist(fullfile(pwd,sprintf('s%d+orig.HEAD',i)),'file')
     unix(cmd);
    end
    
    cmd = sprintf('3dretroicor -prefix pcs%d ',i);
    if ~isempty(resp)
        cmd = sprintf('%s -resp %s -respphase respphase_s%d.1D',cmd,resp{i},i);
    end
    
    if ~isempty(ecg)
        data = load(ecg{i});
      %  threshold = min(data)+(max(data)-min(data))*3/4;
      threshold = 500;  % Assuming the peak height has been normalized to 1000.
      
        fprintf('Scan %d: Threshold %4.3f\n',i,threshold);
        cmd = sprintf('%s -card %s -threshold %f -cardphase cardphase_s%d.1D',cmd,ecg{i},threshold,i);
    end
    cmd = sprintf('%s s%d+orig',cmd,i);
    unix(cmd);
    
  %  cmd = sprintf('3dTshift -rlt+ -prefix s%d_tsh s%d+orig',i,i);
    cmd = sprintf('3dTshift -rlt+ -tzero 0 -prefix pcs%d_tsh pcs%d+orig',i,i);  % aligned to the beginning of each TR 04/21/2010, X.Z.
    %cmd = sprintf('3dcopy s%d+orig s%d_tsh',i,i); % no correction
    if ~exist(fullfile(pwd,sprintf('pcs%d_tsh+orig.HEAD',i)),'file')
      unix(cmd);
    end
    
end    
    
if ~exist('i_volreg','var')
    [err,info] = BrikInfo(name_volreg);
    i_volreg = info.DATASET_RANK(2)-1;
end


% volume register
for i=1:nf
   cmd = sprintf('3dvolreg -base %s''[%d]'' -prefix pcs%d_vr -1Dfile pcmotion_s%d.1D pcs%d_tsh+orig',name_volreg,i_volreg,i,i,i);
   if ~exist(fullfile(pwd,sprintf('pcs%d_vr+orig.HEAD',i)),'file')
    unix(cmd);
   end
   
end

% detrend 
for i=1:nf
   cmd = sprintf('3dDetrend -prefix pcs%d_vr_dtr -polort 2 pcs%d_vr+orig',i,i);
   if ~exist(fullfile(pwd,sprintf('pcs%d_vr_dtr+orig.HEAD',i)),'file')
    unix(cmd);
   end
end


% normalized signal intensity
clipv = zeros(1,nf);
for i=1:nf
    
      cmd =  sprintf('3dTstat -prefix pcs%d_mean pcs%d_vr+orig',i,i);
    if ~exist(fullfile(pwd,sprintf('pcs%d_mean+orig.HEAD',i)),'file')
     unix(cmd); 
    end
   
        cmd = sprintf('3dClipLevel -mfrac 0.55 pcs%d_vr+orig',i);
        [s,clip] = unix(cmd);
        clipv(i)=str2double(clip);
        disp(clip);
        
    if ~exist(fullfile(pwd,sprintf('pcs%d_dtr_norm+orig.HEAD',i)),'file')

        cmd = sprintf('3dcalc -a pcs%d_vr_dtr+orig -b pcs%d_mean+orig -expr ''a/b*step(b-%f)'' -prefix pcs%d_dtr_norm',i,i,clipv(i),i);
        unix(cmd);
    end
   
end

%generate brainmask.
if ~exist(fullfile(pwd,'brainmask+orig.HEAD'),'file')
dsets = '';
expr = sprintf('step(a-%f)',clipv(1));
for i=1:nf
  dsets = sprintf('%s -%c pcs%i_mean+orig',dsets,'a'+i-1,i);
  if i>1
    expr = sprintf('%s*step(%c-%f)',expr,'a'+i-1,clipv(i));
  end
end
cmd = sprintf('3dcalc %s -expr ''%s'' -prefix brainmask',dsets,expr);
unix(cmd);
end

% align anatomy to functional
if ~no_anat && ~exist(fullfile(pwd,'anat_al+orig.HEAD'),'file')
cmd = sprintf('align_epi_anat.py -epi %s -anat anat+orig -epi_base %d -partial_axial',name_volreg,i_volreg);
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
    
    
    
    
