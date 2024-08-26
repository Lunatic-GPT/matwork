a=zeros(64,1);

a([19,34,46])=1;
fa=fft(a);
fm=fftmat(64);

s=randperm(64);
s=s(1:20);

psf=zeros(64);
for i=1:64
    for j=1:64
      ej=zeros(64,1);
      ej(j)=1;
      ei=zeros(1,64);
      ei(i)=1;
      psf(i,j)=ei*fm(s,:)'*fm(s,:)*ej;
    end
end

psf2=psf;
for i=1:64
    psf2(i,i)=0;
end

lambda=0.5;
fftm=fftmat(64);
[m,c1,c2]=cs_recon(fa(s),fftm(s,:),eye(64),lambda);

f2=zeros(1,64);
f2(s)=fa(s);
m2=ifft(f2);

figure;plot(real(m));

figure;plot(real(m2));

fprintf('side lobe max = %f\n',max(abs(psf2(:)))/abs(psf(1,1)));
fprintf('lambda=%f, cost1 = %f, cost2 = %f\n',lambda,c1,c2);
save('simulate_cs_1d.mat');



