function reconTSE_UIH(dname,motionPar,ds)

[n,suf]=strtok(dname,'_');

fname=name4pat(fullfile(dname,['*',suf,'.raw']),1);
if isempty(fname) disp('no file found');return; end

xmlname=name4pat(fullfile(dname,['*',suf,'.prot']),1);
xml=parseXML(xmlname);
FOVro=readPar_uih(xml,'FOVro'); %FOV along [ro, pe, par]
FOVpe=readPar_uih(xml,'FOVpe');

Thickness=readPar_uih(xml,'Thickness');
OverSamplingPE=readPar_uih(xml,'OverSamplingPE');
OverSamplingSPE=readPar_uih(xml,'OverSamplingSPE');
SlicePerSlab=readPar_uih(xml,'SlicePerSlab');
FOV=[FOVro,FOVpe*(1+OverSamplingPE/100),Thickness*SlicePerSlab*(1+OverSamplingSPE/100)];
center_pos=readFOVCenter_uih(xml); %[x,y,z] in dicom convention;

iscan=strtok(dname,'_');
%nro=readPar_uih(xml,'ReadResolution');
%PhaseResolution=readPar_uih(xml,'PhaseResolution');

orient=readOrient_uih(xml);
an=readInplaneRotAngle_uih(xml);

if ~any(orient~=[-1,0,0]) && an==0 %sag
    ro_pe_par=[3,2,-1];
    rotmat=[0,0,1;0,1,0;-1,0,0]';

elseif ~any(orient~=[-1,0,0]) && an==90 %sag
    ro_pe_par=[2,-3,-1];
    rotmat=[0,1,0;0,0,-1;-1,0,0]';
elseif ~any(orient~=[0,-1,0]) && an==0 %cor
    ro_pe_par=[-3,1,-2];
     rotmat=[0,0,-1;1,0,0;0,-1,0]';
elseif ~any(orient~=[0,-1,0]) && an==90 %cor
    ro_pe_par=[1,3,-2];
     rotmat=[1,0,0;0,0,1;0,-1,0]';
elseif ~any(orient~=[0,0,-1]) && an==0 %tra
    ro_pe_par=[-1,2,-3];
     rotmat=[-1,0,0;0,1,0;0,0,-1]';
elseif ~any(orient~=[0,0,-1]) && an==90 %tra
    ro_pe_par=[2,1,-3];
     rotmat=[0,1,0;1,0,0;0,0,-1]';
else
    error('unknown orientation');
end

center_pos=center_pos(abs(ro_pe_par)); %to ro, pe, par

%npe=nro*PhaseResolution/100*FOVpe/FOVro;
%PartialPE=readPar_uih(xml,'PartialPE'); %10 short TE; 11 long TE

[data,~,~,~,~,~,~,shotIndex]= Read_UIH_Raw_v5_7_additional_shot_index(fname);
data=squeeze(data);
data=permute(data,[3,1,2,4]);%change to [ro,pe,par,ch]

shotIndex=shotIndex+1;
shotIndex(squeeze(sos(data(1,:,:,:),4))==0)=0;
voxsize=FOV./size(data(:,:,:,1));


fd=ifft1c(data,1);

clear data;

fd=fd(end/4+1:end-end/4,:,:,:);


method=6;
%2:

if method==1
    %1: motion correction -> nufft' to cartesian grid-> FFT to space -> SAKE -> ESPIRiT; image blur, too many
    %   non-zero entries after nufft; no help even if set positions not in the original sampling mask to zero.  (test results in E:/Projects/MC_TSE/Test_reconTSE_UIH/8_fse_mx_3d_fatnav2_ACS5_reconTSE_UIH_Method1).
    data2=fft1c(fd,1);
    img=recon_method1(data2,motionPar,shotIndex,FOV,ro_pe_par,ds);
elseif method==2
    % SAKE -> ESPIRiT -> motion correction -> nufft
    [img,img_nc,img_c]=recon_method2(fd,motionPar,shotIndex,FOV,ro_pe_par,ds);     

    save(sprintf('%s_recon_Method2.mat',iscan),'img','img_c','img_nc');

    save2nii(img,voxsize,center_pos,rotmat,[iscan,'_espirit']);

    save2nii(img_nc,voxsize,center_pos,rotmat,[iscan,'_espirit_sos']);


    tmp=permute(img_c,abs(ro_pe_par)); %change to [ro,pe,par]

    for i=1:3
      if ro_pe_par(i)<0
        tmp=flip(tmp,i);
      end
    end

    save2nii(tmp,voxsize,center_pos,[iscan,'_espirit_mc_sos']);
    
elseif method==4
% do nufft_espirit at each step of k_ro 
    img=recon_method4(fd,motionPar,shotIndex,FOV,ro_pe_par,ds);
    
save2nii(img,voxsize,center_pos,rotmat,[dname,'_method4']);

elseif method==5
    % motion correction of the fully sampled center lines -> NUFFT3D' to cartesian grid-> SAKE
    % -> ESPIRiT 3D; To be tested at 3T, artifacts at 7 T.
elseif method==6
  % do nufft_spirit at each step of k_ro 
 img=recon_method6(fd,motionPar,shotIndex,FOV,ro_pe_par,ds);

end

function save2nii(img,voxsize,center_pos,rotmat,prefix)

    o.voxsize=voxsize;
    o.center=center_pos; %coord in dicom convention
    o.rotmat=rotmat;
  
    o.sz=size(img);
    o=center2pos_o(o,o.center);

    mat2niigz(img,[],[prefix,'.nii.gz'],true,o);
  
function [img,img_nc,img_c]=recon_method2(fd,motionPar,shotIndex,FOV,ro_pe_par,ds)

%% for test
%load 8_fse_mx_3d_fatnav2_ACS5_fd170.mat
%fd=shiftdim(fd_170,-1);
%%

mask=squeeze(sos(sos(fd,1),4))>1;
[rows,cols]=fullySampledRegion(mask);

sz=size(fd);
img=zeros(sz(1:3));
resESPIRiT=zeros([sz(1:3),2]);
reskESPIRiT=zeros(sz);
weights=zeros([sz(1:3),2]);
%ltic=tic;
parpool(8);
parfor i=1:sz(1)
    dtmp=squeeze(fd(i,:,:,:));
    calib=squeeze(dtmp(rows,cols,:));
    [img(i,:,:),resESPIRiT(i,:,:,:),weights(i,:,:,:),reskESPIRiT(i,:,:,:)]=do_SAKE_ESPIRiT(dtmp,calib);
 %   time_left(i,sz(1),toc(ltic));
end
save recon_method2_reskESPIRiT reskESPIRiT -v7.3;
%%
%  load tmp_testreconTSE.mat
%  motionPar=[];
%  shotIndex=ones(352,270);
[k,data2]=get_k_data_moco_shotIndex(fft1c(reskESPIRiT,1),motionPar,shotIndex,FOV,ro_pe_par,ds);
% k(:,1)=0;


imSize=size(reskESPIRiT(:,:,:,1));
imSize=imSize(abs(ro_pe_par));

im=zeros([imSize,size(data2,2)]);

len=ceil(size(k,1)/10);

for i=1:10
    if i<10
        ind=(i-1)*len+1:i*len;
    else
        ind=(i-1)*len+1:size(k,1);
    end

    nufft=NUFFT3D(k(ind,:),1,[0,0,0],imSize,1,1);
    im=im+nufft'*data2(ind,:);
end



clear data2;
img_c=sos(im,4);
img_nc=sos(ifft1c(ifft1c(reskESPIRiT,2),3),4);

%%

%save 10_recon_images_reg_mid img_c img_nc;

%%




function img=recon_method1(data,motionPar,shotIndex,FOV,ro_pe_par,ds)


mask=squeeze(sos(sos(data,1),4))>1;
[rows,cols]=fullySampledRegion(mask);

[k,data2]=get_k_data_moco_shotIndex(data,motionPar,shotIndex,FOV,ro_pe_par,ds);
sel=sos(data2,2)>0;
k2=k(sel,:);
data2=data2(sel,:);

imSize=size(data(:,:,:,1));
imSize=imSize(abs(ro_pe_par));

nufft=NUFFT3D(k2,1,[0,0,0],imSize,1,1);

im=nufft'*data2;

fim=fft1c(fft1c(im,abs(ro_pe_par(2))),abs(ro_pe_par(3)));

fim=permute(fim,[abs(ro_pe_par),4]); %change to [ro,pe,par,ch];

%rows=size(fim,2)/2-11:size(fim,2)/2+12;
%cols=size(fim,3)/2-11:size(fim,3)/2+12;

%%

sz=size(fd);
img=zeros(sz(1:3));
for i=1:size(fd,1)

    dtmp=squeeze(fim(i,:,:,:));
    calib=squeeze(dtmp(i,rows,cols,:));

    %res=do_SAKE_ESPIRiT(dtmp,calib,ksize,wnthresh,eigThresh_im,nIterCG );
    img(i,:,:)=do_SAKE_ESPIRiT(dtmp,calib);

end


%obsolete 2023/7/12
function m=get_shotindex(d)
% E:/Projects/MC_TSE/Test_reconTSE_UIH/8_fse_mx_3d_fatnav2_ACS5_reconTSE_UIH_Method1).

%the first dimension is shot

sz=size(d);

m=zeros(sz(2:end));
for i=1:sz(1)
    tmp=squeeze(d(i,:,:,:,:));

    if any(vec(m>0 & tmp>0))
        warning('data acquired in multiple shots');
    end
    m(tmp>0)=i;

end







