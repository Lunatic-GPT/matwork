

load('mask_uniform');
m=shiftdim(m,-1);
m=repmat(m,[64,1,1,1]);

run_cs_wavelet_mask2(0.01,0.001,0.001,m,'uniform');




