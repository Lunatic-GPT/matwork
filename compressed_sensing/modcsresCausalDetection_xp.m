function xhat = modcsresCausalDetection_xp(y,mask2,XFM,xfmWeight)

disp('Reconstruting using Modified-CS-Residual')

sz=size(y);
if ndims(mask2)<ndims(y)   % assume 2 d images
    sztmp=ones(1,ndims(y));
    sztmp(end)=sz(end);
    mask2=repmat(mask2,sztmp);
end
y=y./max(y(:));
xhat=zeros(size(y));
res_waveletmodcsres=zeros(size(y));
for seq=1:sz(end)
    
    if seq==1
        
        tmp=sum(y,4);
        ntmp=sum(mask2,4);
        ntmp(ntmp==0)=1;
        y(:,:,:,1)=tmp./ntmp;
        
         FT=p2DFTxp(ones(sz(1:end-1)));
        disp('Reconstructing first frame using inverse DFT:')
        
        xhat_first = FT'*y(:,:,:,1);
                
        wave_first=XFM*xhat_first;
        
        xhat(:,:,:,seq) = xhat_first;
        
    else
        fprintf('Recon image %d\n',seq);
        FT=p2DFTxp(mask2(:,:,:,seq));
        
        mask=mask2(:,:,:,seq);
        
        res0temp=XFM*xhat(:,:,:,seq-1);       
  
        tmp0=conj(res0temp(:)).*res0temp(:);
        tmp=sort(tmp0(:),'descend');
        sumtmp=sum(tmp);
        tmpc=cumsum(tmp);
        
        ind=find(tmpc>=sumtmp*0.99);
        thr_modcsres=sqrt(tmp(ind(1)));
        
        efrac=sum(tmp0(abs(res0temp(:))>thr_modcsres))/sum(tmp);
      
        fprintf('energy fraction is %f; support fraction is %f\n',efrac,length(find(abs(res0temp)>thr_modcsres))/length(tmp));
        
        detectset=find(abs(res0temp)>thr_modcsres);        %known part of set T

       
        res0=zeros(size(mask));  
     
        res0(detectset)=res0temp(detectset);     %Estimation on T
   
        masknz=ones(size(res0));
   
        masknz(detectset)=0;                     %masknz sets all indices in T to be 0, o.w. 1

        yres=y(:,:,:,seq)-FT*xhat_first;        %residual measurements is obtained by substracting the accurate estimation of the first frame, because we are using fullsampling for the first frame
   
          param = init_csparam;
          param.TV=TVOPxp;
          param.TVWeight=0;
          param.FT = FT;
          param.XFM =XFM;
          param.data=yres;
          param.xfmWeight = xfmWeight;
        
          res = XFM*(FT'*yres);                    %Initial guess of the residual signal
          
         
         %%%%Do Modified-CS-residual
         param.masknz=masknz;  %temporary disabled to compare; this is necessary.
          for n=1:2
            
            res = fnlCg_xp(res,param);            %Modified-CS-residual reconstruction
        %%{
          	im_res=param.XFM'*(res+wave_first);
            
            
            figure(101);
            ns=size(im_res,3);
             for isl=1:ns
               subplot(2,ns,isl); imshow(angle(im_res(:,:,isl)'),[]);
          	   subplot(2,ns,isl+ns); imshow(abs(im_res(:,:,isl)'),[]), drawnow;
             end
            %}
          end 
          
          res_waveletmodcsres(:,:,:,seq)=wave_first+res;  % reconstructed results=residual+estimation
          
          xhat(:,:,:,seq)=XFM'*res_waveletmodcsres(:,:,:,seq);  %reconstructed image
    
   end
                  
end
   