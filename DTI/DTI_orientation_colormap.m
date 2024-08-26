function DTI_orientation_colormap(d,nrow,ncol,dc)
% 

d=d./sos(d,4);
do_flip=0;
if do_flip
d=flip(permute(d,[2,1,3,4]),1);
dc=flip(permute(dc,[2,1,3,4]),1);

end

d=ri(d);
d=abs(d);

FA=calc_FA(dc);

d=uint8(d.*FA*255);

sz=size(d);

im=reshape(d,[sz(1),sz(2),ncol,nrow,sz(4)]);
im=permute(im,[1,4,2,3,5]);
im=reshape(im,[sz(1)*nrow,sz(2)*ncol,sz(4)]);

figure;imshow(im);