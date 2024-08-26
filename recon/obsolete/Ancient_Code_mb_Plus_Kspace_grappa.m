%% 
elseif pat==2
    
    Ks=-3:3;
    Ksy=[-2,0,2];
 %   iph=ind_max(2)-lc/2:2:ind_max(2)+lc/2-1;
 %   kernel1=sliceGRAPPAKernel(dsb2(iread,iph,:,:,1),Ks,Ksy);
    iph=ind_max(2)-lc/2:ind_max(2)+lc/2-1;
    kernel=sliceGRAPPAKernel(dsb2(iread,iph,:,:,1),Ks,Ksy);
    
  %  kernel=(kernel1+kernel2)/2;
  %   iph=ind_max(2)-lc/2+1:ind_max(2)+lc/2;
   %  kernel_new=sliceGRAPPAKernel(dsb2(iread,iph,:,:,1),Ks,Ksy);
   
     
   end

   
   
   %%
   %% test GRAPPA = 2 recon
if pat==2
    kernelX=[-2,-1,0,1,2];
    kernelY=[-3,-1,1,3];
    
    
    for iv=1:2
        for sl=1:nsl
            kernel_grp=myGRAPPA_2data_xp([],dsb2(iread,iph,:,sl,1),1:length(iph),1:length(iph),kernelY,kernelX,0,false);
            
           if lstart>1
                lindex2=lstart-1:pat:size(dmb2,2);
           else
                lindex2=lstart+1:pat:size(dmb2,2);
           end
            tmp=myGRAPPA_2data_xp(res(:,:,:,sl,iv),[],lindex,[],kernelY,kernelX,kernel_grp,false);
            res(:,lindex2,:,sl,iv)=tmp(:,lindex2,:);
            
        end
    end
    
end

%%

%%
%{ 

%GRAPPA is not good for this situation because the neibouring k-space line
is critical for successful unaliasing of slices which is not acquired with
GRAPPA.



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
