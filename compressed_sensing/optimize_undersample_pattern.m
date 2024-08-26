a=fftmat(64);
a=circshift(a,[32,32]);
nt=100;

 ref=load('refg_0_10_30_6cy_TR2.0');
%m=mask_cs1d(fid_prefix);

m=cstable_1dgauss(64,4,1,nt,1);
%m=cstable_gems_mTR(64,nt,4,1,6,'');
m=cstable_gems_mTR(64,nt,4,1,0,'');

%a=rand(n,nt);
%m=a>0.75;

m=ones(64,1,nt);
nnz=sum(m(:,1,1));
X=zeros(nt*nnz,64*2);

for i=0:nt-1
%    m2=repmat(m(:,j,i+1),[1,64]);

    m2=repmat(m(:,1,i+1),[1,64]);
    
    X(i*nnz+1:(i+1)*nnz,:)=reshape(a(m2>0),[nnz,64]);
end
b=X'*X;
disp(trace(inv(b)));


    
    
    