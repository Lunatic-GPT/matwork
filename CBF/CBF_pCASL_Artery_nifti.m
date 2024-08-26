function CBF_pCASL_Artery_nifti(fname,tlab,tpostlab,pc,T1a,eff)
% res=CBF_pCASL_Artery(ds[,tlab,tpostlab,pc,T1a,eff])
% fname: the nifti file from dcm2niix
% tlab: labeling time
% tpostlab: postlab delay
% pc: tissue-blood partition coefficient
% T1a: artery T1 in s
% eff: labelling efficiency
% result in ml/100g/min



if ~exist('pc','var')
    pc=0.9; %unit ml/g
end

if ~exist('eff','var')|| isempty(eff)
    eff=0.68;
end

if ~exist('T1a','var')|| isempty(T1a)
    T1a=1/0.67; %unit ml/g
end

if ~exist('tpostlab','var')|| isempty(tpostlab)
    tpostlab=1; %adFree[2]/1000000
end

if ~exist('tlab','var')|| isempty(tlab)
    tlab=18.5*82/1000; % adFree[3]*18.5/1000
end

nii=load_untouch_niigz(fname);
sc=double(nii.img(:,:,:,2:2:end));
sl=double(nii.img(:,:,:,1:2:end-1));

ds=1-sl./sc;

nii.img=CBF_pCASL_Artery(ds,tlab,tpostlab,pc,T1a,eff);

save_untouch_niigz(nii,['CBF_',strtok(filename(fname),'.')]);

if size(nii.img,4)>1
    nii.img=mean(nii.img,4);
    save_untouch_niigz(nii,['CBFmean_',strtok(filename(fname),'.')]);
end









