load('mask_uniform_nc6');
m=shiftdim(m,-1);
m=repmat(m,[64,1,1,1]);

nois=0.01;
tmpname = sprintf('synthesize_kdata_mask2_nois%4.3f_bold%4.3f.mat',nois,0.02);
tmp=load(tmpname);
z2=tmp.z3;

z2=m.*z2;
z2(:,:,:,1)=tmp.z3(:,:,:,1);
XFM=Wavelet_ISU([64,64]);
xfmWeight=0.1;
xhat=modcsresCausalDetection_xp(z2,m,XFM,xfmWeight);

prefix = sprintf('modCS_uniform_nc6_xfmW%4.3f',xfmWeight);
 
write_afni(abs(xhat),prefix);
write_afni(angle(xhat),[prefix,'_ph']);
        
        