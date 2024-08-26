nois_array=[0.02,0.04,0.10];

for j=1
    
 nois=nois_array(j);
 tmpname = sprintf('synthesize_kdata_mask2_nih_nois%4.3f_bold%4.3f.mat',nois,0.015);
 tmp=load(tmpname);
 z2=tmp.z3;


 pat={'uniform','uniform_nc6','gauss'};
 for i=1

   m=get_mask(pat{i});
        
   
         for k=1:size(z2,3)            
          fprintf('Slice %d\n',k);
          tmp=z2(:,:,k,:).*m(:,:,k,:);
          tmp=permute(tmp,[2,1,3,4]);
          tmp=circshift(tmp,[size(tmp,1)/2,size(tmp,2)/2,0,0]);
          tmp=squeeze(tmp);
          recon(:,:,k,:)=cart_ktFOCUSS_KTFOCUSS(tmp);        
         end
 
     recon=circshift(recon,[size(tmp,1)/2,size(tmp,2)/2,0,0]); 
     recon=permute(recon,[2,1,3,4]);
   %write_afni(abs(recon),sprintf('ktFOCUSS_nih_nois%4.3f_%s',nois,pat{i}));

end

end

