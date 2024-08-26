%%
load('mask_uniform');
m=shiftdim(m,-1);
m=repmat(m,[64,1,1,1]);
run_cs_ft_mask2(0.01,0.002,0.002,m,'uniform');