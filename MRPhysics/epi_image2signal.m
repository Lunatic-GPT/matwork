function epi_image2signal(fname,brik)
opt.Frames = brik;
d = BrikLoad(fname,opt);
sz = size(d);
fd = zeros(size(d));
for z=1:size(d,3)
    
    fd(:,:,z) = fft2(fftshift(d(:,:,z)));
    
    fd(:,:,z) = fftshift(fd(:,:,z));
    
    sig = zeros(1,sz(1)*sz(2));
    
    for j=1:sz(2)
        if mod(j,2) == 0
            sig((j-1)*sz(1)+1:j*sz(1)) = fd(:,j,z);
        else
            sig((j-1)*sz(1)+1:j*sz(1)) = flipdim(fd(:,j,z),1);
        end
    end
figure;
plot(1:sz(1)*sz(2),abs(sig));

end




