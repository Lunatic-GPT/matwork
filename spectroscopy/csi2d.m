function csi2d(fid)

a=read_fid(fid);
seqcon=readPar(fid,'seqcon');
sw = readPar(fid,'sw');


lsfid=14;
if strcmp(seqcon(2:end-1),'ncssn');
    nv=readPar(fid,'nv');
    nv2=readPar(fid,'nv2');
    np=readPar(fid,'np');
    a=reshape(a,[np/2,nv,nv2]);
end

tof=0;
np2=np/2-lsfid;
freq=((1:np2)-np2/2-1)*sw/np2+tof;
freq=freq/400+4.675;

fa=fft(a,[],2);
fa=fft(fa,[],3);
ph=repmat(angle(fa(lsfid+1,:,:)),[size(fa,1),1,1]);
t = 0:size(a,1)-lsfid-1;
t = t./sw;
wt=exp(-(t/0.05).^2);
wt = repmat(wt',[1,nv,nv2]);
fa=fa.*exp(-1i*ph);

fa=ifft(fa(lsfid+1:end,:,:).*wt,[],1);
fa=fftshift(fa,1);
%fa=circshift(fa,size(a)/2);
fa = fftshift(fa,2);
fa = fftshift(fa,3);


figure;

n1=size(a,2);
n2=size(a,3);
for i=n1/2+1:n1
    for j=1:n2
         disp([i,j]);
      subplot(n2,n1/2,(j-1)*n1/2+i-n1/2);
      plot(freq,abs(fa(:,i,j)));
      ind=(freq>1&freq<6);
      ylim([min(abs(fa(ind,i,j))),max(abs(fa(ind,i,j)))]);
      xlim([1,6]);
      hold on;
     % plot(imag(fa(:,i,j)),'r');
      
      box off;
      xlabel('');
      ylabel('');
      set(gca,'XTick',1:4,'YTick',[]);
      set(gca,'Xdir','reverse');
    end
end
set(gcf,'Position',[370,12,1542,976]);

      