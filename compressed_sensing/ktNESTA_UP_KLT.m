function xt=ktNESTA_UP_KLT(z,m,pred)
% adapt from ktNESTA_UP, not yet implemented.
% xt=ktNESTA_UP(z,m,pred)
% z and m are a matrix with dimension of nx*ny*1*nt
% m is a mtrix with dimension of nx*ny*1*nt
 scl=max(abs(z(:)));
 
 z=z/scl;
 
if ~exist('pred','var')  
    pred=0;
    zpred=0;
end

if pred
    z2=z;
    z2(m==0)=0;
    z2=sum(z2,4);
    m2=sum(m,4);
    
    z2=z2./m2;
    z2(isnan(z2))=0;
    zpred=repmat(z2,[1,1,1,size(z,4)]);
    z=z-zpred;
else
    zpred=0;
end

    

A=@(x) kt_mat(x,m);
At=@(x) kt_matt(x,m);

 Ac = @(x) counter( @(x) kt_mat(x,m), x);
 Atc= @(x) counter( @(x) kt_matt(x,m), x);
    

 
sigma=0.001;
n=numel(m);

Lambda=sigma*sqrt(2*log(n));
  La = my_normest( A, At, n)^2;
% the factor was 0.1
muf=0.1*sigma/La;
opts = [];
opts.MaxIntIter = 6;
opts.TOlVar = 1e-6;
opts.verbose = 50;
opts.maxiter = 3000;

opts.U=speye(n);
opts.Ut=speye(n);
opts.stoptest = 1;  
opts.errFcn=@(x) errFcn(x,z,m);


z2=z(m>0);

xf=NESTA_UP(Ac,Atc,z2(:),Lambda,La,muf,opts);



xf=reshape(xf,size(m));
%xf=xf+repmat(ifft2c(zpred(:,:,:,1)),[1,1,1,size(zpred,4)]);


xt=fft1c(xf,4)*scl+scl*repmat(ifft2c(zpred(:,:,:,1)),[1,1,1,size(zpred,4)]);


%save ktNESTA_up;


function out=errFcn(x,z,m)
      kt=kt_mat(x,m);
      out=norm(kt-z(m>0))/sum(m(:));
      
function out=kt_mat(b,m)
        % b in xf space
        % out in kt space
     
        sz=size(m);
        b=reshape(b,sz);
  
        c= fft1c(b,4);
        
        fc=fft2c(c);
   
        out=fc(m>0); 
        out=out(:);
        
function out=kt_matt(b,m)
     % b in kt space
     % out in xf space
    
        d=zeros(size(m));
        
        d(m>0)=b;
        out=ifft1c(ifft2c(d),4);
        out=out(:);
        
        