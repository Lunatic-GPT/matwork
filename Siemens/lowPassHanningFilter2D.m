function im_lowpass=lowPassHanningFilter2D(im,n)


fim=fft2c(im);

sz=size(im);
if length(sz)==2
    sz(3)=1;
end

im_lowpass=0*im;

    filter=0*fim;
    
    [xx,yy]=meshgrid(0:n-1,0:n-1);
    
    filter_mat=0.25*(1-cos(2*pi*xx/(n-1))).*(1-cos(2*pi*yy/(n-1)));
    
    
    filter(sz(1)/2-n/2+1:sz(1)/2+n/2,sz(2)/2-n/2+1:sz(2)/2+n/2,:,:)=repmat(filter_mat,[1,1,sz(3),size(im,4)]);
    
    
    
    

    im_lowpass=ifft2c(fim.*filter);
