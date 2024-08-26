function [im,data2]=recon_grappa2D(data,kindex,kmax,data_cal,kindex_cal)
%[im,kdata]=recon_grappa2D(data,kindex,kmax,[data_cal,kindex_cal])
%data should be in the format of nro*npe_us*nslice*ncoil*extra
% kindex should be of length npe_us; 1 based.
% kmax, the total number of k-indice after recon.
% data_cal: in the format of nro_cal*npe_us_cal*nslice*ncoil*extra; only the first
% of extra dimension will be used.
% kindex_cal should be of length npe_us_cal; 1 based. kindex and kindex_cal
% im: the output dim: nro*kmax*nslice*extra


 sz=size(data);
 sz(2)=kmax;
 data2=single(zeros(sz));
 data2(:,kindex,:,:,:)=data;
 kindex2=sort(kindex);
 acs=acsLines(kindex2);
 
 if exist('data_cal','var')
 
     kindex_cal2=sort(kindex_cal);
     acs_cal=acsLines(kindex_cal2);
 end
 
 
 kernelX=[-2,-1,0,1,2];
 if abs(kindex2(2)-kindex2(1))==3
     kernelY=[-5,-2,1,4;-4,-1,2,5];
 elseif abs(kindex2(2)-kindex2(1))==2
    kernelY=[-3,-1,1,3];
 elseif abs(kindex2(2)-kindex2(1))==4
     kernelY=[-6,-2,2,-6;-7,-3,1,5;-5,-1,3,7];
 end
 
 %%
 isampl=kindex2(1):diff(kindex2(1:2)):kmax;
 
 %%
 
 for isl=1:size(data,3)
     tic;
     for ik=1:size(kernelY,1)
         
         if exist('data_cal','var')
             sz=size(squeeze(data_cal(:,:,1,:)));
             sz(2)=max(kindex_cal2);
             
             tmp=zeros(sz);
             tmp(:,kindex_cal2,:)=data_cal(:,:,isl,:,1);
             [tmp_empty,coef]=myGRAPPA_xp(tmp,kindex_cal2,acs_cal,kernelY(ik,:),kernelX,0,false);
         end
         
         for i5=1:size(data,5)
             if exist('data_cal','var') 
                 Rk = myGRAPPA_xp(squeeze(data2(:,:,isl,:,i5)),kindex2,acs,kernelY(ik,:),kernelX,coef);
             else
                 Rk = myGRAPPA_xp(squeeze(data2(:,:,isl,:,i5)),kindex2,acs,kernelY(ik,:),kernelX,0);
             end
             

             ineg_kern=kernelY(ik,kernelY(ik,:)<0);
             ineg_kern=sort(ineg_kern);
             
             isampl_index=isampl-ineg_kern(end);
             isampl_index(isampl_index>kmax)=[];
             data2(:,isampl_index,isl,:,i5)=Rk(:,isampl_index,:,:,:);
             
         end
     end
     time_left(isl,size(data,3),toc);
 end
 
 sz=size(data2);
 sz=[sz(1:3),sz(5:end)];
 
 im=zeros(sz,'single');
 for i=1:size(data2,3)

         for k=1:size(data2,5)
          im(:,:,i,k)=sos(ifft2c(data2(:,:,i,:,k)),4);        
              
         end
 end
 
     
function a=acsLines(sk)


a=[];
for i=1:length(sk)-1
    if abs(sk(i)-sk(i+1))==1
        a(end+1)=sk(i);
    elseif i>1&&abs(sk(i)-sk(i-1))==1
        a(end+1)=sk(i);
    end
end
