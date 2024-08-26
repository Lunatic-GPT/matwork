function [kdata,m2,frequency]=kdata_FLASH2DCS_Zspectrum(a)
%m2=kdata_FLASH2DCS(a)


d=readb_fid(a);


pe=readbPar(fullfile(a,'acqp'),'ACQ_spatial_phase_1');
NR = readbPar(fullfile(a,'method'),'PVM_NRepetitions');
%riset=D(4);

%rampt = D(5);%readbPar(fullfile(a,'method'),'PVM_RampTime');

%askip = round((riset+rampt)*sw);
mtrx=readbPar(fullfile(a,'method'),'PVM_Matrix');
pe=round(pe*mtrx(1)/2)+mtrx(1)/2+1;
nav=d(1,:,:,:);
data=d(end-mtrx(1)+1:end,:,:,:);

data2=phaseCorr(data,nav);

if mod(NR,length(pe)/mtrx(2))~=0
   warning('mod(NR,length(pe)/mtrx(2))~=0');
end

f_pTR = readbPar(fullfile(a,'method'),'f_pTR');

reduction=mtrx(2)/f_pTR;
data2=reshape(data2,[mtrx(1),f_pTR,NR]);

pe=reshape(pe,[mtrx(2)/reduction,length(pe)*reduction/mtrx(2)]);
%pe=pe';
of1=readbPar(fullfile(a,'method'),'PVM_SPackArrPhase1Offset');
of0=readbPar(fullfile(a,'method'),'PVM_SPackArrReadOffset');
fov=readbPar(fullfile(a,'method'),'PVM_Fov');

sz=size(data);

d2=zeros([mtrx',NR]);
m2=zeros(mtrx(1),mtrx(2),1,NR);

%dataf=readb_fid('15');
%dataf=dataf(end-63:end,:,:);
%dataf=fov_shift_kspace(dataf,[0,-of1],fov);

for j=1:NR
    j2=mod(j-1,size(pe,2))+1;
  
  d2(:,pe(:,j2),j)=data2(:,:,j);
  %d2(:,pe(:,j2),j)=dataf(:,pe(:,j2),j);  
  m2(:,pe(:,j2),j)=1;
end


frf=readbPar(fullfile(a,'method'),'freqlist');

z=reshape(d2,[mtrx',1,length(frf),NR/length(frf)]);
z=mean(z,5);
kdata=fov_shift_kspace(z,[0,-of1],fov);  % comment out if use a fully sampled data.

m2=reshape(m2,[mtrx',1,NR]);
m2=m2(:,:,:,1:length(frf));

%%
frequency=readbPar(fullfile(a,'method'),'freqlist');


%frf2=frf(abs(frf)<10000);
%m2=m2(:,:,:,abs(frf)<10000);
%z2=z(:,:,:,abs(frf)<10000);
if nargout==0
  save(sprintf('%s_RawData',a),'kdata','m2','frequency');
  if reduction==1
    fz=fft2c(z);
    write_afni(abs(fz),sprintf('%s_recon_phasecorr',a));
    return;
  end


end


function d2=phaseCorr(d2,nav)

  ph=angle(nav);
  
  ph=repmat(ph,[size(d2,1),1,1,1]);
  
 d2=d2.*exp(-1i*ph);

