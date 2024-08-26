function reconLS_FLASH(a,format)
%reconLS_FLASH(a,format)

if ~exist('format','var')
    format='s';
end
d=readb_fid(a);

sw=readbPar(fullfile(a,'method'),'PVM_EffSWh');
D=readbPar(fullfile(a,'acqp'),'D');
pe=readbPar(fullfile(a,'acqp'),'ACQ_spatial_phase_1');
NR = readbPar(fullfile(a,'method'),'PVM_NRepetitions');
kzero = find(pe==0);
riset=D(4);

rampt = D(5);%readbPar(fullfile(a,'method'),'PVM_RampTime');

askip = round((riset+rampt)*sw);
mtrx=readbPar(fullfile(a,'method'),'PVM_Matrix');
pe=round(pe*mtrx(1)/2)+mtrx(1)/2+1;
nav=d(askip+1:askip+mtrx(1),:,:,:);
data=d(end-mtrx(1)+1:end,:,:,:);
data2=data;
if mod(NR,length(pe)/mtrx(2))~=0
   error('mod(NR,length(pe)/mtrx(2))~=0');
end

iNR = length(pe)/mtrx(2);
data=reshape(data,[mtrx(1),iNR,mtrx(2),NR/iNR]);
data=permute(data,[1,3,2,4]);

nav=reshape(nav,[mtrx(1),iNR,mtrx(2),NR/iNR]);
nav=permute(nav,[1,3,2,4]);

tmp=reshape(data,[mtrx',1,NR]);
%write_afni(angle(tmp),[a,'_kspace_ph']);

%data=phaseCorr(data,nav);

nav=reshape(nav,[mtrx',1,NR]);
%write_afni(angle(nav),[a,'_navkspace_ph']);


tmp=reshape(data,[mtrx',1,NR]);
%write_afni(angle(tmp),[a,'_kspace_ph_corr']);

%write_afni(abs(nav),[a,'_nav']);
%write_afni(angle(nav),[a,'_nav_ph']);

%fnav=fft1c(nav,1);
%write_afni(abs(fnav),[a,'_fnav']);
%write_afni(angle(fnav),[a,'_fnav_ph']);

pe=reshape(pe,[iNR,mtrx(2)]);
pe=pe';

d2=zeros(size(data));
for j=1:NR/iNR
 for i=1:iNR
  d2(:,pe(:,i),i,j)=data(:,:,i,j);
 end
end

d2=reshape(d2,[mtrx',1,NR]);


of1=readbPar(fullfile(a,'method'),'PVM_SPackArrPhase1Offset');
of0=readbPar(fullfile(a,'method'),'PVM_SPackArrReadOffset');
fov=readbPar(fullfile(a,'method'),'PVM_Fov');
d2=fov_shift_kspace(d2,[0,-of1],fov);

%write_afni(angle(d2),[a,'_kspace_ph_corr']);
img=fft2c(d2);

img=flipdim(img,1);
img=flipdim(img,2);

if format=='a'
  %  write_afni(abs(img),[a,'_recon_nav']);
    write_afni(angle(img),[a,'_recon_ph']);
elseif format=='s'
    writesdt4(abs(img),[a,'_recon']);
else
    error('unknown format');
end

figure;
plot(abs(nav(:,kzero(1))));
hold on;
y=flipdim(data2(:,kzero(1)),1);

plot(abs(y),'r');

legend('Nav Echo','k=0 data');

function d2=phaseCorr(d2,nav)

for i=1:size(nav,2)
  ph=phase(squeeze(nav(32,i,:)));
  
  ph2=lowpass(ph,50);
%ph2=ph;
  ph2=shiftdim(ph2,-2);
  ph2=repmat(ph2,[size(d2,1),1,1]);
  mag=abs(nav);
  mag=abs(nav);

 d2(:,i,:)=d2(:,i,:).*exp(-1i*ph2*8/5.5);
end




  
