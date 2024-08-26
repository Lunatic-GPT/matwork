function recon_ESPIRiT3D_T2_core(fname,dfile,fmap_ind,prefix_map,ine,file_res0,debug)

% fname='meas_MID85_tse_vfl_pss_FatNav_FID23810.mat';
% dfile = 'Motion_MID85.1D';
% fov = [210.0,165.7,97.6];  % FOV before accounting for the oversampling
% factor
% lMatrix = [512,405,244,256];
% 8/25/2019: remove fov and lMatrix parameters;


load(fname);

Data = fft1c(Data,1);   % to k space;
if ~exist('debug','var')
    debug=false;
end

prefix=strtok2(fname,'.');
pro=fullfile(prefix,[prefix,'.pro']);
fov(1)=readsPar(pro,'asSlice[0].dReadoutFOV');
fov(3)=readsPar(pro,'asSlice[0].dThickness');
fov(2)=readsPar(pro,'asSlice[0].dPhaseFOV');


lBaseResolution = readsPar(pro,'lBaseResolution');
lPhaseEncodingLines =readsPar(pro,'lPhaseEncodingLines');
lPartitions = readsPar(pro,'lPartitions');
lImagesPerSlab = readsPar(pro,'lImagesPerSlab');

lMatrix=[lBaseResolution,lPhaseEncodingLines,lImagesPerSlab,lPartitions];

if debug
   [Data,Line,Partition,lMatrix] = crop_data(Data,Line,Partition,lMatrix); 
end

pnav_sag = readsPar(pro,'sCuboid.sPosition.dSag');
pnav_cor= readsPar(pro,'sCuboid.sPosition.dCor');
pnav_tra=readsPar(pro,'sCuboid.sPosition.dTra');
if isempty(pnav_sag)
    pnav_sag=0;
end
if isempty(pnav_cor)
    pnav_cor=0;
end
if isempty(pnav_tra)
    pnav_tra=0;
end

p_sag = readsPar(pro,'asSlice[0].sPosition.dSag');
p_cor= readsPar(pro,'asSlice[0].sPosition.dCor');
p_tra=readsPar(pro,'asSlice[0].sPosition.dTra');


if isempty(p_sag)
    p_sag=0;
end
if isempty(p_cor)
    p_cor=0;
end
if isempty(p_tra)
    p_tra=0;
end

nro=size(Data,1);
Nc=32;
%

fov_os=[fov(1:2),fov(3)*lMatrix(4)/lMatrix(3)];  % fov_os after oversampling along Partition direction
ro_pe_par=[1,-2,-3];
ds = [p_sag,p_cor,p_tra]-[pnav_sag,pnav_cor,pnav_tra];
[k2,Data]=get_k_data(Data,Line,Partition,dfile,Nc,fov_os,lMatrix([2,4]),ro_pe_par,ds);



imSize=load(fmap_ind,'imSize');
imSize=imSize.imSize;


nufft=NUFFT3DwMap(k2,1,[0,0,0],imSize,prefix_map,fmap_ind,ine);

dfile2=filename(dfile);
if ~isempty(dfile2)
prefix_save=unique_name([prefix,'_',strtok(dfile2,'.'),'_',num2str(ine)]);
else
prefix_save=unique_name([prefix,'_',num2str(ine)]);
    
end

if ~exist('file_res0','var') || isempty(file_res0)
    
    disp('Calculating initial image');
    res=nufft'*Data;
    disp('Initial image calculation done');
    
    save([prefix_save,'_iter0.mat'],'res', '-v7.3');
    res0=res;
    
else
    res0=load(file_res0);
    res0=res0.res;
end

XOP = Wavelet('Daubechies',4,6);
disp('Start CGL1ESPIRiT');
tcg=tic;
nIterCG = 5;
res = cgL1ESPIRiT_nomap(Data, res0, nufft, nIterCG,XOP,0,0,1);
disp('CGL1ESPIRiT finished');
toc(tcg);

if exist('file_res0','var') && ~isempty(file_res0)
    prefix_save=file_res0(1:end-4);
end
save([prefix_save,'_iter',num2str(nIterCG),'.mat'], 'res','-v7.3');



function [Data,Line,Partition,l] = crop_data(Data,Line,Partition,lMatrix) 

    
    l=[64,64,lMatrix(3)/lMatrix(4)*32,32];
    Line=int16(Line);
    Partition=int16(Partition);
    c=floor(lMatrix/2);
    sel=(Line-c(2))>=-l(2)/2 & (Line-c(2))<l(2)/2 & (Partition - c(4))>=-l(4)/2 & (Partition - c(4))<l(4)/2;
    ind=size(Data,1)/2-l(1)/2+1:size(Data,1)/2+l(1)/2;
    Data=Data(ind,sel);
    
    Line=Line(sel)-c(2)+l(2)/2;
    Partition=Partition(sel)-c(4)+l(4)/2;