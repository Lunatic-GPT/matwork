function recon_ESPIRiT3D_T2(fname,dfile,nmapjobs,file_res0,ine,debug)

% fname='meas_MID85_tse_vfl_pss_FatNav_FID23810.mat';
% dfile = 'Motion_MID85.1D';
% fov = [210.0,165.7,97.6];  % FOV before accounting for the oversampling
% factor
% lMatrix = [512,525,244,256];
% 8/25/2019: remove fov and lMatrix parameters;
% ine: the dimension not to carry out espirit;
% i.e. if ine ==1, then espirit will be carried out along pe and par
%      if ine == 3, then espirit will be carried out along ro and pe; the
%      par dimension will be zero filled for partial Fourier acquirsition before convert to image domain.
%

if ~exist('debug','var')
    debug=false;
end

prefix=strtok2(fname,'.');


fov(1)=readsPar(fullfile(prefix,[prefix,'.pro']),'asSlice[0].dReadoutFOV');
fov(3)=readsPar(fullfile(prefix,[prefix,'.pro']),'asSlice[0].dThickness');
fov(2)=readsPar(fullfile(prefix,[prefix,'.pro']),'asSlice[0].dPhaseFOV');


lBaseResolution = readsPar(fullfile(prefix,[prefix,'.pro']),'lBaseResolution');
lPhaseEncodingLines =readsPar(fullfile(prefix,[prefix,'.pro']),'lPhaseEncodingLines');
lPartitions = readsPar(fullfile(prefix,[prefix,'.pro']),'lPartitions');
lImagesPerSlab = readsPar(fullfile(prefix,[prefix,'.pro']),'lImagesPerSlab');

lMatrix=[lBaseResolution,lPhaseEncodingLines,lImagesPerSlab,lPartitions];

pnav_sag = readsPar(fullfile(fname(1:end-4),[prefix,'.pro']),'sCuboid.sPosition.dSag');
pnav_cor= readsPar(fullfile(fname(1:end-4),[prefix,'.pro']),'sCuboid.sPosition.dCor');
pnav_tra=readsPar(fullfile(fname(1:end-4),[prefix,'.pro']),'sCuboid.sPosition.dTra');
if isempty(pnav_sag)
    pnav_sag=0;
end
if isempty(pnav_cor)
    pnav_cor=0;
end
if isempty(pnav_tra)
    pnav_tra=0;
end

p_sag = readsPar(fullfile(fname(1:end-4),[prefix,'.pro']),'asSlice[0].sPosition.dSag');
p_cor= readsPar(fullfile(fname(1:end-4),[prefix,'.pro']),'asSlice[0].sPosition.dCor');
p_tra=readsPar(fullfile(fname(1:end-4),[prefix,'.pro']),'asSlice[0].sPosition.dTra');


if isempty(p_sag)
    p_sag=0;
end
if isempty(p_cor)
    p_cor=0;
end
if isempty(p_tra)
    p_tra=0;
end

load(fname);

Data = fft1c(Data,1);   % to k space;
if debug
    [Data,Line,Partition,lMatrix] = crop_data(Data,Line,Partition,lMatrix);
end

nro=size(Data,1);
Nc=32;

%


disp('calculate maps');
fov_os=[fov(1:2),fov(3)*lMatrix(4)/lMatrix(3)];  % fov_os after oversampling along Partition direction
ro_pe_par=[1,-2,-3];
ds = [p_sag,p_cor,p_tra]-[pnav_sag,pnav_cor,pnav_tra];
imSize=lMatrix([1,2,4]);
do_nufft3D_calib=1;
if do_nufft3D_calib
    [k2,Data]=get_k_data(Data,Line,Partition,dfile,Nc,fov_os,imSize(2:3),ro_pe_par,ds);
else
    [Data2,mask]=sort_kData(Data,Line,Partition,Nc,imSize);  %only for debug;
    Data3=crop(Data2,[512,24,24,32]);
end

%prep_bsub_ESPIRiT3D_maps(fname,dfile,fov_os,lMatrix([2,4]),nmapjobs,ncomp);
%% old prep_bsub_ESPIRiT3D_maps

lRefLinesPE=readsPar(fullfile(fname(1:end-4),[prefix,'.pro']),'sPat.lRefLinesPE'); % fix for now



nmaps = 2;

ind=linspace(0,imSize(ine),nmapjobs+1);
ind=round(ind);
ind1=ind(1:end-1)+1;
ind2=ind(2:end);

mid=strtok_no(fname,'_',2);
if ~isempty(dfile)
    out_prefix=[mid,'_',strtok(filename(dfile),'.'),'_ine',num2str(ine)];
else
    out_prefix=[mid,'_ine',num2str(ine)];
end

fmap_ind = sprintf('%s_Params_%d.mat',out_prefix,ine);
save(fmap_ind,'ind1','ind2','nro','Nc','nmaps','imSize');


if ispc
    prefix_map=prep_bsub_ESPIRiT3D_maps(k2,Data,Nc,imSize,lRefLinesPE,out_prefix,nmapjobs,nmaps,ine);
    recon_ESPIRiT3D_T2_core(fname,dfile,fmap_ind,prefix_map,ine,file_res0,debug);
else
    prefix_map=prep_bsub_ESPIRiT3D_maps_sbatch(k2,Data,Nc,imSize,lRefLinesPE,out_prefix,nmapjobs,nmaps,ine);
    wait_for_files(prefix_map,ind1,ind2);
    
    script_name=sprintf('recon_%s.sh',out_prefix);
    
    func=sprintf('recon_ESPIRiT3D_T2_core(''%s'',''%s'',''%s'',''%s'',%d,''%s'',%d)',fname,dfile,fmap_ind,prefix_map,ine,file_res0,debug);
    make_matlab_shellscript(script_name,func,400,48);
    cmd=sprintf('source %s',script_name);
    unix(cmd);
    
end

function [Data,Line,Partition,l] = crop_data(Data,Line,Partition,lMatrix)

Line=int16(Line);
Partition=int16(Partition);
l=[64,64,lMatrix(3)/lMatrix(4)*32,32];

c=floor(lMatrix/2);
sel=(Line-c(2))>=-l(2)/2 & (Line-c(2))<l(2)/2 & (Partition - c(4))>=-l(4)/2 & (Partition - c(4))<l(4)/2;
ind=size(Data,1)/2-l(1)/2+1:size(Data,1)/2+l(1)/2;
Data=Data(ind,sel);

Line=Line(sel)-c(2)+l(2)/2;
Partition=Partition(sel)-c(4)+l(4)/2;


function [Data3,mask]=sort_kData(Data2,Line,Partition,Nc,imSize)

p=Line(1:Nc:end)-min(Line);
s=Partition(1:Nc:end)-min(Partition);

Data2=reshape(Data2,[size(Data2,1),Nc,size(Data2,2)/Nc]);

nskip=max(find(diff(double(p))==0))+1;% last index for phase correction scan
if isempty(nskip) % no phase correction scan; only skip the GRAPPA noise scan
    nskip=1;
end

Data3 = zeros([imSize(1),imSize(2),imSize(3),Nc],'single');
mask=zeros([imSize(2),imSize(3)]);
%%

p2 = p(nskip+1:end);  % this needs to be changed.
s2=s(nskip+1:end);

for i=1:length(s2)
    Data3(:,p2(i)+1,s2(i)+1,:)=Data2(:,:,nskip+i);
    mask(p2(i)+1,s2(i)+1) = 1;
end

function fname2=save_calib_data(prefix,k2,Data,Nc,imSize,lRefLinesPE,ine)

sel=abs(k2(:,2))<=54/imSize(2)/2 & abs(k2(:,3))<=54/imSize(3)/2;
Data = reshape(Data,[length(Data(:))/Nc,Nc]);

k2_center=k2(sel,:);
Data=Data(sel,:);

imSize_lowres = [54,54,54];%[length(pe_center),length(ro_center),length(par_center)];
imSize_lowres(ine)=imSize(ine);  % no change in the readout direction.


nufft=NUFFT3D(k2_center.*repmat(floor(imSize/2)./floor(imSize_lowres/2),[size(k2_center,1),1]),1,[0,0,0],imSize_lowres,1,1);

im = nufft'*Data;
ift=[1,2,3];
ift(ine)=[];
kData=fft1c(fft1c(im,ift(1)),ift(2));   % this matches with d2

acs_line=floor(size(kData,2)/2)-lRefLinesPE/2+1:floor(size(kData,2)/2)+lRefLinesPE/2;

if ine==1
    acs_par=size(kData,3)/2-15:size(kData,3)/2+16;
    kCalib=squeeze(kData(:,acs_line,acs_par,:));
else
    acs_ro=size(kData,1)/2-15:size(kData,1)/2+16;
    kCalib=squeeze(kData(acs_ro,acs_line,:,:));
end
fname2 = sprintf('%s_kCalib.mat',prefix);
save(fname2, 'kCalib','imSize');


function prefix2=prep_bsub_ESPIRiT3D_maps_sbatch(k2,Data,Nc,imSize,lRefLinesPE,prefix,njobs,nmaps,ine)

fname=save_calib_data(prefix,k2,Data,Nc,imSize,lRefLinesPE,ine);

ind=linspace(0,imSize(ine),njobs+1);
ind=round(ind);

prefix2=sprintf('%s_maps',strtok(fname,'.'));

for i=1:length(ind)-1
    
    fmapname=sprintf('%s_%d_%d.mat',prefix2,ind(i)+1,ind(i+1));  %do not regenerate
    if exist(fmapname,'file')
        continue;
    end
    
    
    func=sprintf('calc_ESPIRiT3D_maps(''%s'',%d:%d,%d,%d)',fname,ind(i)+1,ind(i+1),nmaps,ine);
    
    fname2=unique_name(sprintf('%s_%d.sh',prefix2,ind(i)+1));
    make_matlab_shellscript(fname2,func,32,24);
    cmd=sprintf('source %s',fname2);
    disp(['submit ',fname2]);
    pause(60);  % pause so that all scripts will be submitted properly, hopefully.
    unix(cmd);
    
end



function prefix2=prep_bsub_ESPIRiT3D_maps(k2,Data,Nc,imSize,lRefLinesPE,prefix,njobs,nmaps,ine)

fname=save_calib_data(prefix,k2,Data,Nc,imSize,lRefLinesPE,ine);

ind=linspace(0,imSize(ine),njobs+1);
ind=round(ind);

prefix2=sprintf('%s_maps',strtok(fname,'.'));
for i=1:length(ind)-1
    
    fmapname=sprintf('%s_%d_%d.mat',prefix2,ind(i)+1,ind(i+1));  %do not regenerate
    if exist(fmapname,'file')
        continue;
    end
    
    calc_ESPIRiT3D_maps(fname,ind(i)+1:ind(i+1),nmaps,ine);%d:%d,%d) logmaps_%d\n',fname2,ind(i)+1,ind(i+1),nmaps,ind(i)+1)
    
end



