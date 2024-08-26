function mat=fft2mat(n,m)

m1=fftmat(n);
m2=fftmat(m);
m2=conj(m2');

mat=m1(:)*m2(:)';
        
   
