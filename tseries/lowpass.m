function out=lowpass(ts,np,dim)
% out=lowpass(ts,np,dim)
% np: number of points to keep nonzero, on each side of the k=0 point.
if ~exist('dim','var')

f=fft(ts);
%f=fftshift(f);

f(np+2:end-np)=0;

out=ifft(f);

else
    
f=fft(ts,[],dim);
%f=fftshift(f);
if dim==1
    f(np+2:end-np,:,:,:)=0;
elseif dim==2
 f(:,np+2:end-np,:,:)=0;
elseif dim==3
 f(:,:,np+2:end-np,:)=0;
elseif dim==4
    f(:,:,:,np+2:end-np)=0;
end
    
out=ifft(f,[],dim);
end
