function [im,data2,coef_all]=recon_grappa2D_Siemens(data,kindex,kmax)
%[im,kdata]=recon_grappa2D(data,kindex,kmax,[data_cal,kindex_cal])
%data should be in the format of nro*npe_us*nslice*ncoil*extra
% both nro and nslice dimensions should be already in the image space.
% kindex should be of length npe_us; 1 based.
% kmax, the total number of k-indice after recon.
% data_cal: in the format of nro_cal*npe_us_cal*nslice*ncoil*extra; only the first
% of extra dimension will be used.
% kindex_cal should be of length npe_us_cal; 1 based. kindex and kindex_cal
% im: the output dim: nro*kmax*nslice*extra
% in Siemens ice: the ro dimension was IFT'ed before grappa.

 sz=size(data);
 sz(2)=kmax;
 data2=single(zeros(sz));
 data2(:,kindex,:,:,:)=data;
 kindex2=sort(kindex);
 acs=acsLines(kindex2);
 
 if exist('data_cal','var')
 
     kindex_cal2=sort(kindex_cal);
     acs_cal=acsLines(kindex_cal2);
     data_cal=ifft1c(data_cal,1);
 end
 
 
kernelX=0;
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
 lx=length(kernelX);
 ly=size(kernelY,2);
 Ncoils=size(data,4);
coef_all = zeros(size(data,1),lx*ly*Ncoils,size(kernelY,1),Ncoils,size(data,3),size(data,5));

    tm= tic;
for iro=1:size(data,1)
   
 for isl=1:size(data,3)
 
       for ik=1:size(kernelY,1)
             
         for i5=1:size(data,5)
             
             data2_tmp=reshape(data2(iro,:,isl,:,i5),[1,size(data2,2),size(data2,4)]);
             
             [Rk,coef] = myGRAPPA_xp(data2_tmp,kindex2,acs,kernelY(ik,:),kernelX,0);
             

             ineg_kern=kernelY(ik,kernelY(ik,:)<0);
             ineg_kern=sort(ineg_kern);
             
             isampl_index=isampl-ineg_kern(end);
             isampl_index(isampl_index>kmax)=[];
             data2(iro,isampl_index,isl,:,i5)=Rk(1,isampl_index,:,:,:);
             coef_all(iro,:,ik,:,isl,i5)=coef;
         end
   
       end
           
 end
   time_left(iro,size(data,1),toc(tm));
end
 
 sz=size(data2);
 sz=[sz(1:3),sz(5:end)];
 
 im=zeros(sz,'single');
 for i=1:size(data2,3)
         for k=1:size(data2,5)
          im(:,:,i,k)=sos(ifft1c(data2(:,:,i,:,k),2),4);        
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
