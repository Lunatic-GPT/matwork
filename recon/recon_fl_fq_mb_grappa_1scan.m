function recon_fl_fq_mb_grappa_1scan(fmb)

%fmb='meas_MID42_fl_fqZ_AF3_FID10324.dat';
%fsb='meas_MID44_fl_fqZ_SB_FID10326.dat';

% GRAPPA recon

[d,lin,par,sl]=readMeasDat(fmb,inf,0,true);

prefix_mb=strtok(fmb,'.');
nsl=readsPar([prefix_mb,'.pro'],'alFree[4]');
nave=readsPar([prefix_mb,'.pro'],'lAverages');

nline=length(lin)/(nsl+nave)/64;

dsb=d(:,1:64*nsl*nline);
dmb=d(:,64*nsl*nline+1:end);
lin_mb=lin(:,64*nsl*nline+1:end);

lin2=lin_mb(1:64*nave:end);



af=readsPar([prefix_mb,'.pro'],'alFree[6]');
nave_sb=1;

nro=size(dmb,1);
dmb=reshape(dmb,[nro,32,2,nave,length(dmb(:))/nro/64/nave]);
dsb=reshape(dsb,[nro,32,2,nave_sb,length(dsb(:))/nro/64/nsl/nave_sb,nsl]);

%%

dmb2=zeros(nro,32,2,max(lin2(:))+1,'single');

dsb2=zeros(nro,32,2,max(lin2(:))+1,nsl,'single');

    for i=1:length(lin2)
        dmb2(:,:,:,lin2(i)+1)=mean(dmb(:,:,:,:,i),4);
        dsb2(:,:,:,lin2(i)+1,:)=mean(dsb(:,:,:,:,i,:),4);
    end

dsb2=interleave2linear(dsb2,5);

fall=zeros(nro,32,2,nline,nsl);
for i=1:nsl
    
f=exp(1i*(0:nline-1)*2*pi/af*(i-1-(nsl-1)/2));
f=reshape(f,[1,1,1,nline]);

fall(:,:,:,:,i)=repmat(f,[nro,32,2,1]);

end

%%

im_sb=ifft1c(dsb2,4);
im_sb=permute(im_sb,[1,4,2,5,3]);
im_sb=squeeze(sos(im_sb,3));
save(['reconSBRef_',prefix_mb,'.mat'],'im_sb');


fdsb2=fft1c(dsb2,1);   % k space
fdsb2=fdsb2.*fall;

fdmb2=fft1c(dmb2,1);  %k space

% determine the partial fourie factor along RO
tmp=sos(fdmb2,2);
[tmp2,ind_max]=max(tmp(:));
ind_max=ind2subb(size(tmp),ind_max);
negphase=ind_max(1)-1;

% GRAPPA
fdsb2=permute(fdsb2,[1,4,5,3,2]);
fdmb2=permute(fdmb2,[1,4,5,3,2]);

iread=floor(negphase/2):floor(negphase/2)*3;
disp('Slice Grappa Calibration');
kernel=sliceGRAPPAKernel(permute(squeeze(fdsb2(iread,end/8*3:end/8*5,:,1,:)),[1,2,4,3]),3,3);

res=zeros(nro,nline,32,nsl,2,'single');
Ks=3;
disp('Slice Grappa Recon');
for i=1:2
res(:,:,:,:,i)=sliceGRAPPAKernel(squeeze(fdmb2(:,:,:,i,:)),Ks,Ks,kernel);
end
save(['kdataSliceGrappa_',prefix_mb],'res','kernel');

disp('Partial FT recon');
%%{
im_mb=0*res(:,Ks+1:end-Ks,:,:,:);
for isl=1:nsl
    
        for j=1:2
            
            tmp=res(:,Ks+1:end-Ks,:,isl,j);
            
     %       tmp=partialFT(tmp(:,1:negphase+nro/2),negphase);
          
             kdata2=cat(1,zeros(nro/2-negphase,nline-2*Ks,32),tmp(1:negphase+nro/2,:,:)*10000);
             
             
             %kdata2=permute(kdata2,[3,1,2]);

             tmp=partialFT_pocs(permute(kdata2,[3,1,2]),20,true);
            
            im_mb(:,:,:,isl,j)=permute(tmp,[2,3,1]);
            
            disp([isl,j]);
        end
    
end
%}
%{
im_mb=0*res;
for isl=1:nsl
    for i=1:32
        for j=1:2
            tmp=permute(res(:,:,i,isl,j),[2,1]);           
            tmp=partialFT(tmp(:,1:negphase+nro/2),negphase);
           %  kdata2=cat(1,zeros(320-negphase,640,32),kdata(1:negphase+320,:,:));
           %  kdata2=permute(kdata2,[3,1,2]);
           %  tmp=partialFT_pocs(tmp(:,1:negphase+nro/2),negphase);
            im_mb(:,:,i,isl,j)=permute(tmp,[2,1,3,4]);
            disp([isl,i,j]);
        end
    end
end

%}
%%

im_mb2=ifft2c(fft2c(im_mb).*permute(conj(fall(:,:,:,Ks+1:end-Ks,:)),[1,4,2,5,3]));

%im_mb2=fftshift(fftshift(im_mb2,1),2); %only needed for Cecil's partialFT
pc_mb=angle(mean(im_mb2(:,:,:,:,1).*conj(im_mb2(:,:,:,:,2)),3))*180/pi;
pc_mb=squeeze(pc_mb);
save(['recon_PFTpocs_',prefix_mb],'im_mb2');

save(['reconPC_PFTpocs_',prefix_mb],'pc_mb');

tmp=squeeze(sos(im_mb2,3));
save(['reconMag_PFTpocs_',prefix_mb],'tmp');


%%
%{ 

%GRAPPA is not good for this situation

fdmb2_shift=fftshift(fdmb2,2);  % the first element is k=0;

dsb3=ifft1c(fdsb2,2);  % image space along the pe; -N/2 to N/2-1
dsb3=fftshift(dsb3,2); % 0 to N-1 (or equivalentl 0 to N/2-1; -N/2 to -1)
dsb3=reshape(dsb3,[nro,nline*nsl,2,32]); 
fdsb3_shift=fft(dsb3,[],2);   % the first element is k=0;


fdmb3_shift=0*fdsb3_shift;

fdmb3_shift(:,1:nsl:end,:,:)=fdmb2_shift;

fdsb3=fftshift(fdsb3_shift,2);
fdmb3=fftshift(fdmb3_shift,2);

%%

 kernelX=[-2,-1,0,1,2];
 if nsl==3
     kernelY=[-5,-2,1,4;-4,-1,2,5];
 elseif nsl==2
    kernelY=[-3,-1,1,3];
 elseif nsl==4
     kernelY=[-6,-2,2,-6;-7,-3,1,5;-5,-1,3,7];
 end
 
 
 for i=1
     
     k_sel=1;
     
     for iknly=1:size(kernelY,1)
         coef=myGRAPPA_2data_xp([],squeeze(fdsb3(:,:,k_sel,:)),1:nsl:nsl*nline,nsl*nline/2-15:nsl*nline/2+16,kernelY(iknly,:),kernelX,0);
         res=myGRAPPA_2data_xp(squeeze(fdsb3(:,:,i,:)),[],1:nsl:nsl*nline,nsl*nline/2-15:nsl*nline/2+16,kernelY(iknly,:),kernelX,coef);
       %  coef=myGRAPPA_2data_xp([],fdsb3,1:nsl:nline,nline/2-15:nline/2+16,kernelY(iknly,:),kernelX,0);
       %  res=myGRAPPA_2data_xp(fdmb3,[],1:nsl:nline,nline/2-15:nline/2+16,kernelY(iknly,:),kernelX,coef);
       
        %  res=myGRAPPA_xp(squeeze(fdsb3(:,:,i,:)),1:nsl:nsl*nline,nsl*nline/2-15:nsl*nline/2+16,kernelY(iknly,:),kernelX,0);
      %     res=myGRAPPA_xp(squeeze(fdsb3),1:nsl:nline,nline/2-15:nline/2+16,kernelY(iknly,:),kernelX,0);
          
           ineg_kern=kernelY(iknly,kernelY(iknly,:)<0);
          ineg_kern=sort(ineg_kern);
          isampl_index=1:nsl:nline*nsl;
          
          isampl_index=isampl_index-ineg_kern(end);
          isampl_index(isampl_index>nline*nsl)=[];
          isampl_index(isampl_index<1)=[];
          
          fdmb3(:,isampl_index,i,:)=res(:,isampl_index,:);
          
          disp([i,iknly]);
     end
     
     
 end
 
 tmp=ifft2c(fdmb3);
 
 figure;imshow(sos(tmp(:,:,1,:),4),[]);
%}
 %%
%

%%


