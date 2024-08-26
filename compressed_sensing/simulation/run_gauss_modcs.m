load('mask_gauss');
m=permute(m,[1,2,4,3]);


nois=0.01;
tmpname = sprintf('synthesize_kdata_mask2_nih_nois%4.3f_bold%4.3f.mat',nois,0.02);
tmp=load(tmpname);
z2=tmp.z3;

z2=m.*z2;
z2(:,:,:,1)=tmp.z3(:,:,:,1);
XFM=Wavelet_ISU([64,64]);
xfmWeight=0.1;
xhat=modcsresCausalDetection_xp(z2,m,XFM,xfmWeight);

prefix = sprintf('modCS_nih_gauss_xfmW%4.3f',xfmWeight);
 
write_afni(abs(xhat),prefix);
write_afni(angle(xhat),[prefix,'_ph']);
        
        