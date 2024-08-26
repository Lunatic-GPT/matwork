%z2=synthesize_kdata(0.005,0.2);

%tmp=load('synthesize_kdata_nois0p005_BOLD0p04');
tmp=load('synthesize_kdata_nois0p005_BOLD0p04_sin');
z2=tmp.z2;
recon_orig=ifft2c(z2);
a=loadtable('petable_R4_N64_nTR240_ns1');
m=zeros(64,64,1,240);
for i=1:240
    m(:,a((i-1)*16+1:i*16)+33,1,i)=1;
end

z2=m.*z2;
recon=zeros(size(z2));
        for i=1:size(z2,3)            
          fprintf('Slice %d\n',i);
         %{
          tmp=z2(:,:,i,:);
          tmp=permute(tmp,[2,1,3,4]);
          tmp=circshift(tmp,[size(tmp,1)/2,size(tmp,2)/2,0,0]);
          tmp=squeeze(tmp);
          recon(:,:,i,:)=cart_ktFOCUSS_KTFOCUSS(tmp);
          recon(:,:,i,:)=circshift(recon(:,:,i,:),[size(z2,1)/2,size(z2,2)/2,0,0]);
          recon(:,:,i,:)=permute(recon(:,:,i,:),[2,1,3,4]);
     
          %}
          %%{
           ftt=FTt;
          [recon(:,:,i,:),d,e]=ktFOCUSS_xp(z2,m,0,1/2,ftt*recon_orig);  
      %}
        end
        recon=recon*64;
  %      write_afni(abs(recon),'test_ktFOCUSS');
  %      write_afni(abs(recon_orig),'recon_orig');
        
 %       save ktFOCUSS_default recon 
  % recon=ifft2c(z2);
        mask=zeros(64,64);
        mask(48:52,27:31)=1;
      %  mask(50,29)=1;
        ts=mean_roi(abs(recon),mask);
        ts_orig=mean_roi(abs(recon_orig),mask);
        
        figure;subplot(2,1,1);plot(ts);
        
        subplot(2,1,2);plot(ts_orig);
        
        figure;subplot(1,2,1);imshow(abs(recon(:,:,1,1)),[]);
        subplot(1,2,2);imshow(abs(recon_orig(:,:,1,1)),[]);