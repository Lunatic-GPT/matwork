N=8;
a=fftmat(N);
sa=mat2supmat(a,N,'r');
sb=mat2supmat(conj(a'),N,'l');
sab=sb*sa;
isab=inv(sab);

dx=Dxmat(N);
dy=dx';

sdx=mat2supmat(dx,N,'r');
sdy=mat2supmat(dy,N,'l');


Lambda=(sdx'*sdx+sdy'*sdy)*isab*ones(N*N,1);

reshape(Lambda,[N,N])
for i=1:size(Lambda,1)
 %   Lambda(i,i)=0;
end

%%
ifft2(ones(N,N))
uker2=reshape(diag(Lambda),[N,N]);
dx'*dx*ifft2(ones(N,N))+ifft2(ones(N,N))*dy*dy'
%

%
rows=16;
cols=16;

uker(1,1) = 4;uker(1,2)=-1;uker(2,1)=-1;uker(rows,1)=-1;uker(1,cols)=-1;
uker = fft2(uker);