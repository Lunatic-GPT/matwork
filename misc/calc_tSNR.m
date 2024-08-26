function calc_tSNR(fname)

nii=load_untouch_niigz(fname);
d=single(nii.img);
rSNR=mean(d,4)./std(d,[],4);
nii.img=rSNR;


prefix=strtok(filename(fname),'.');

save_untouch_niigz(nii,['tSNR_',prefix]);







