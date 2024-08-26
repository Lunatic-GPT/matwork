function mask_PAinPC_SiemensDCM(fpat,scanIDs,wm_mask,alpha_2tail,prefix)
% mask_PAinPC(fpat,scanIDs,wm_mask,interp_factor)
% 4/29/2021:
% calcu a PA mask based on a series of phase and mag images
% tSNR calculated to determine threshold for roi determination.
% par should contain a parameter named alpha_2tail 
% changed name from mask_PAinPC to mask_PAinPC_SiemensDCM

warning off;
nscans=length(scanIDs);

for i=1:nscans
    dname=sprintf(fpat,scanIDs(i)); 
    data(:,:,:,:,i)=ri(dname,1);    
end

data(:,:,:,1,:)=data(:,:,:,1,:)/4095*2*pi;
voxsize = dcmDimCenter(dname);

sz=size(data);
mag=reshape(data(:,:,:,2,:),[sz(1:3),5]);
ph=reshape(data(:,:,:,1,:),[sz(1:3),5]);

roi=mask_PAinPC(mag,ph,wm_mask,alpha_2tail,voxsize,prefix);

save('prefix','roi');
