function xt2=ktSENSE_NESTA_UP(z,m,sigma,pred)
% xt=ktNESTA_UP(z,m[,pred])
% z and m are matrices with dimension of ny*nt*nchan and ny*nt
% calculate prediction; default false.

 scl=1;%max(abs(z(:)));
 
 z=z/scl;
 
if ~exist('pred','var')  
    pred=0;
   
end

nchan=size(z,3);
nt=size(z,2);

    z2=z;
    m2=repmat(m,[1,1,nchan]);
    z2(m2==0)=0;
    z2=sum(z2,2);
    m2=sum(m2,2);
    
    z2=z2./m2;
    z2(isnan(z2))=0;
    z0=z2;
if pred
    zpred=repmat(z0,[1,size(z,2),1]);
    z=z-zpred;
end
   
img0=ifft1c(z0,1);
sens=img0./repmat(sos(img0),[1,1,32]);


A=@(x) kt_mat(x,m,sens);
At=@(x) kt_matt(x,m,sens);

 Ac = @(x) counter( @(x) kt_mat(x,m,sens), x);
 Atc= @(x) counter( @(x) kt_matt(x,m,sens), x);
     
n=numel(m);

Lambda=sigma*sqrt(2*log(n));
  La = my_normest( A, At, n)^2;
% the factor was 0.1
muf=0.01*sigma/La;
opts = [];
opts.MaxIntIter = 6;
opts.TOlVar = 1e-5;
opts.verbose = 50;
opts.maxiter = 1000;

opts.U=speye(n);
opts.Ut=speye(n);
opts.stoptest = 1;  
opts.errFcn=@(x) errFcn(x,z,m,img0);


z2=z(repmat(m,[1,1,nchan])>0);

xf=NESTA_UP(Ac,Atc,z2(:),Lambda,La,muf,opts);



xf=reshape(xf,size(m));
%xf=xf+repmat(ifft2c(zpred(:,:,:,1)),[1,1,1,size(zpred,4)]);

xt=fft1c(xf,2);


xt2=zeros(size(z));

for i=1:nt
    if pred
     xt2(:,i,:)=repmat(xt(:,i),[1,1,nchan]).*sens+img0;     
    else
      xt2(:,i,:)=repmat(xt(:,i),[1,1,nchan]).*sens;
    end
end

xt2=sos(xt2)*scl;
end

%save ktNESTA_up;
function out=errFcn(x,z,m,img0)
      kt=kt_mat(x,m,img0);
      nchan=size(img0,3);
      m=repmat(m,[1,1,nchan]);
      out=norm(kt-z(m>0))/sum(m(:));      
end

function out=kt_mat(b,m,img0)
        % b in xf space; ny*nt
        % m; ny*nt
        % out in kt space; (ny/r)*nt*nchan
        % img0 in xt space; averaged k-space data; ny*1*nchan
        sz=size(m);
        b=reshape(b,sz);
       
        c= fft1c(b,2);
      %  c=fft1c(c,2);
      %  c=b;
        nchan=size(img0,3);
        nt=size(b,2);
        ny=size(img0,1);
        d=zeros(ny,nt,nchan);
        
        for i=1:nt
            for j=1:nchan
              d(:,i,j)=img0(:,j).*c(:,i);
            %  d(:,i,j)=c(:,i);
            end
        end
        
        fd=fft1c(d,1);
      % fd=d;
        m=repmat(m,[1,1,nchan]);
        out=fd(m>0); 
        out=out(:);
        
end

function out=kt_matt(b,m,img0)
     % b in kt space; (ny/r)*nt*nchan
     % out in xf space; ny*nt
     % m; ny*nt
     % img0 in xt space; averaged k-space data; ny*1*nchan
     
     sz=size(m);
     nchan=size(img0,3);
     ny=size(img0,1);
     nt=size(m,2);
     
     d=zeros([sz,nchan]);
        m=repmat(m,[1,1,nchan]);
        d(m>0)=b;
    %    fd=ifft(ifft(d,[],1),[],2);  %ny*nt*nchan
     
      fd=ifft1c(d,1);
    
   %  fd=ifft1c(fd,2);
     
        out=zeros(ny,nt);
        
        for i=1:nt    
             out(:,i)=sum(conj(img0).*fd(:,i,:),3);
         %      out(:,i)=sum(fd(:,i,:),3);
        end
        
       
    out=ifft1c(out,2);
        out=out(:);
        
        
end
        
        
        