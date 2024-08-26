load('mask_gauss');
m=permute(m,[1,2,4,3]);

run_cs_wavelet_mask2(0.01,0.001,0.001,m,'gauss');




