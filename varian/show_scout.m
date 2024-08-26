function show_scout(fid_prefix,i_first)

if ~exist('i_first','var')
    i_first=1;
end
    
a=read_fid([fid_prefix,'.fid/fid']);
a=dcCorr(a);

z_shft = fftshift(a);
fz = fft2(z_shft);
fz2 = fftshift(fz);

ni = size(a,3);
im = abs(fz2);


figure;

npanel = ni-i_first+1;
if npanel>3
  npanel =3;
end

for i=1:npanel
 subplot(1,npanel,i);
 imagesc(im(:,:,i_first+i-1));
 colormap(gray);
end



function z = dcCorr(z_tmp)
        
 NUM_DC_SUM = 8;
 ave1= mean(z_tmp(1:NUM_DC_SUM/2,:,:,:),1);
 ave2= mean(z_tmp(end-NUM_DC_SUM/2+1:end,:,:,:),1);

 ave=(ave1+ave2)/2;
        
 ave3 = mean(ave,2);       
 z=z_tmp-repmat(ave3,[size(z_tmp,1),size(z_tmp,2),1,1]);
        
