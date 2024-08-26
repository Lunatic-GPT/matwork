function [im,fd]=recon_grappa2D_Lustig(data,kindex,kmax)
% [im,fd]=recon_grappa2D_Lustig(data,kindex,kmax)
%data should be in the format of nro*npe_us*nslice*ncoil
% kindex should be of length npe_us; 1 based.
% kmax, the number of k-indice without under-sampling.

 sz=size(data);
 sz(2)=kmax;
 data2=zeros(sz);
 data2(:,kindex,:,:,:)=data;
 
 mask=zeros(sz(1:2));
 mask(:,kindex)=1;
 coils=size(data,4);
  [tmp,dcomp]=getCalibSize(mask);
  kSize=[5,5];
  for i=1:size(data,3)
    datacomp = squeeze(data2(:,:,i,:)).*repmat(dcomp,[1,1,coils]);
scale_fctr = norm(datacomp(:))/sqrt(coils)/20;
datatmp = squeeze(data2(:,:,i,:))/scale_fctr;

 acs=acsLines(kindex);
 n=size(data2,1);
 kCalib=datatmp(n/2-19:n/2+20,acs,:);
 
  res = GRAPPA(datatmp,kCalib,kSize,0.01);
  end
  
 fd=fft2c(res);
 
  %im=squeeze(sqrt(mean(abs(fd(:,:,:,17:32,:)).^2,4)));
 im=squeeze(sqrt(mean(abs(fd).^2,4)));

 
 
     
function a=acsLines(sk)


a=[];
for i=1:length(sk)-1
    if abs(sk(i)-sk(i+1))==1
        a(end+1)=sk(i);
    elseif i>1&&abs(sk(i)-sk(i-1))==1
        a(end+1)=sk(i);
    end
end
