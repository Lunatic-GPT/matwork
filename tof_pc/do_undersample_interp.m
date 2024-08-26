function res=do_undersample_interp(img,undersample,interp)
% undersample is a fraction, e.g. 0.25 for 4 times undersampling

f=fft2c(img);

ld=round(size(f(:,:,1))*undersample);

sout=round(size(f(:,:,1))*undersample*interp);

if size(f,4)>1    
   f=crop(f,ld(1),ld(2),size(f,3),size(f,4));
elseif size(f,3)>1
   f=crop(f,ld(1),ld(2),size(f,3));
else
   f=crop(f,ld(1),ld(2));
end

f=zpad(f,sout(1),sout(2),size(f,3),size(f,4));

res=ifft2c(f);
