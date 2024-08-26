dname = '11_t2_mx_sag_A10_iso0.8mm_FatNav';
motionPar = 'reg_res/reg_11_t2_mx_sag_A10_iso0.8mm_FatNav.txt';
ds =[0,0,0];

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


%nro=readPar_uih(xml,'ReadResolution');
%PhaseResolution=readPar_uih(xml,'PhaseResolution');

orient=readOrient_uih(xml);
an=readInplaneRotAngle_uih(xml);

if ~any(orient~=[-1,0,0]) && an==0 %sag
    ro_pe_par=[3,2,-1];
    rotmat=[0,0,1;0,1,0;-1,0,0]';

elseif ~any(orient~=[-1,0,0]) && an==90 %sag
    ro_pe_par=[2,-3,-1];
elseif ~any(orient~=[0,-1,0]) && an==0 %cor
    ro_pe_par=[-3,1,-2];
elseif ~any(orient~=[0,-1,0]) && an==90 %cor
    ro_pe_par=[1,3,-2];
elseif ~any(orient~=[0,0,-1]) && an==0 %tra
    ro_pe_par=[-1,2,-3];
elseif ~any(orient~=[0,0,-1]) && an==90 %tra
    ro_pe_par=[2,1,-3];
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



fd=ifft1c(data,1);

clear data;

fd=fd(end/4+1:end-end/4,:,:,:);

voxsize=FOV./size(fd(:,:,:,1));
method=4;
%%
ifd=fft1c(fd,1);
img=recon_method4(ifd,motionPar,shotIndex,FOV,ro_pe_par,ds);


