function fm_showphase(fid_prefix)

a=read_fid([fid_prefix,'.fid/fid']);

bs=(mean(a(1:10,:),1)+mean(a(end-10:end,:),1))/2;
a=squeeze(a);
a=a-repmat(bs,[size(a,1),1]);
 
gs = readPar(fid_prefix,'gf');
gsf = readPar(fid_prefix,'gfs');
lp = readPar(fid_prefix,'lp');
rp = readPar(fid_prefix,'rp');
sw=readPar(fid_prefix,'sw');
lro=readPar(fid_prefix,'lro');
dt=1/sw;

a=squeeze(a);
s=zeros(256,12);

for i=1:size(a,2)
  s(:,i)=wft(a(:,i),gs,gsf,dt,rp,lp,sw);
end

figure;
x = linspace(-lro/2,lro/2,size(s,1));
for i=1:6
  subplot(3,2,i);
  da = angle(s(:,i*2-1))-angle(s(:,i*2));
  da=da*180/pi;
  plot(x,da,'k-');
  ylim([-360,360]);
  xlabel('x (cm)');
end

figure;
x = linspace(-lro/2,lro/2,size(s,1));
for i=1:6
  subplot(3,2,i);
  
% plot(1:size(s,1),angle(s(:,i*2-1))*180/pi);
  plot(x,abs(s(:,i*2-1)),'k-');
  hold on;
 % plot(x,real(s(:,i*2)),'r-');
 % ylim([-180,180]);
  xlabel('x (cm)');
  
set(gca,'XDir','reverse');
end


function s=wft(fid,gs,gsf,dt,rp,lp,sw)


t=(0:length(fid)-1)*dt;
t=t';
f= fid.*exp(-((t-gsf)/gs).^2);
%f=fftshift(f);

s=fft(f,256);
s=fftshift(s);
s=sp_phasecorr(s,rp,lp,sw);


