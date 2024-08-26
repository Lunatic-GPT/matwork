function d2=Oversampling_rm_FT(d)

ftmpData=ifft1c(d,1);

len=length(ftmpData);

l4=round(len/4);

d2=ftmpData(l4+1:l4+len/2,:);

d2=fft1c(d2,1);



