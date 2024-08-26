load('mask_gauss');
m=permute(m,[1,2,4,3]);
run_cs_ft_mask2(0.01,0.002,0.002,m,'gauss_test2');




