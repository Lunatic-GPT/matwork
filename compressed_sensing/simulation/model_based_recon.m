load('c:/users/xiaopeng/Desktop/cs_simulation/synthesize_kdata_mask2_nih_nois0.010_bold0.015.mat');
load('c:/users/xiaopeng/Desktop/cs_simulation/mask_gauss.mat');
ref=load('refg_0_10_30_6cy_TR2.0');

ref2=[ones(1,length(ref));ref(:)'];
xfm=tTransform(ref2);
TVWeight=0;
xfmWeight=0;
nodisplay=false;
m3=permute(m,[1,2,4,3]);
im_res=run_cs(z3,z_true,[],[],m3,xfm,TVWeight,xfmWeight,nodisplay);