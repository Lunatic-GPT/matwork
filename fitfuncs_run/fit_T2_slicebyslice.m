function fit_T2_slicebyslice(slc)

d=ri('Ice_images_zShimmed.mat');

TE =  [0,15,25,35,53.9];
T2=T2map_data(d(:,:,slc,:),TE);

str=sprintf('_%d',slc);
save(['T2p',str], 'T2');

