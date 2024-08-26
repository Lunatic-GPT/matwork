function res = times(a,b)

if isa(a,'fft2o') == 0
    error('In  A.*B only A can be fft2o operator');
end

if a.adjoint
    %res = ifft(b,[],1);
    %res = ifft(res,[],2);
    res=ifft2c(b);
else
    %res = fft(b,[],1);
    %res = fft(res,[],2);
    res=fft2c(b);
end


