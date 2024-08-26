nois_array=[0.02,0.04,0.10];

bold = 0.015;
for j=1
    
 nois=nois_array(j);
 tmpname = sprintf('synthesize_kdata_mask2_nih_nois%4.3f_bold%4.3f.mat',nois,bold);
 tmp=load(tmpname);
 z2=tmp.z3;


 pat={'uniform','uniform_nc6','gauss'};
 for i=3
        
   m=get_mask(pat{i});
        
 z2=m.*z2;
 z2(:,:,:,1)=tmp.z3(:,:,:,1);
 XFM=Wavelet_ISU([64,64]);
 xfmWarray=[0.05,0.1,0.05];
 xfmWeight=xfmWarray(i);
 xhat=modcsresCausalDetection_xp(z2,m,XFM,xfmWeight);

 prefix = sprintf('modCS_nih_%s_nois%f_bold%4.3f_xfmW%4.3f',pat{i},nois,bold,xfmWeight);
 
 write_afni(abs(xhat),prefix);
 write_afni(angle(xhat),[prefix,'_ph']);
        
end

end