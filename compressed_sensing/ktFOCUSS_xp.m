function [res,d,e,f]=ktFOCUSS_xp(v,m,lambda,p,rho)
% not working yet need to figure out how to do matrix inversion
% no need for matrix inversion.  using a conjugate gradient method to
% optimize the cost function.

if ~exist('p','var')
    p=1/2;
end

sz=size(m);
FT=p2DFTxp(m,sz(1:2),1,2);
ftt=FTt;
v=v.*m;
if ~exist('rho','var')
 
   in=sum(v,4);
   
   nm=sum(m,4);
   m2=(nm==size(m,4));
   m2=repmat(m2,[1,1,1,size(m,4)]);
   
   v2=zeros(size(v));
   v2(m2)=v(m2);
   rho=FT'*(ftt*v2);

end

%rho0=rho;
%rho=zeros(size(rho));
%v=v-ft*(ftt'*rho0);

    W=abs(rho).^p; 
    q=rho./W;
for i=1:20
    d=(v-FT*(ftt'*rho)).^2;
    d=sum(abs(d(:)));
    e=rho./W;
    e=sum(abs(e(:).^2));
    
    f=sum(abs(rho(:)));
    rms=sqrt(d/sum(m(:)))./max(abs(v(:)));
    fprintf('Cost = %f; rms = %f; data = %f; reg = %f; sparsity cost =%f. \n',d+e*lambda,rms,d,e,f);
    
    W=abs(rho).^p;   
    rho=W./q;
    
    
end

res=ftt'*(rho);

return;
cov=get_cov(res);
ftt=klt(cov);

%mres=mean(res,4);
%mres=repmat(mres,[1,1,1,size(rho,4)]);
mres=zeros(size(res));
v=v-ft*mres;
rho=ftt'*(res-mres);    
return;
for i=1:10
   
    W=abs(rho).^p;    
    
    d=(v-ft*(ftt*rho)).^2;
    d=sum(abs(d(:)));
    e=lambda*rho./W;
    e=sum(abs(e(:)));
  
    rms=sqrt(d)/sum(m(:))./max(abs(v(:)));
    % the cost function is 
    fprintf('Cost = %f; rms = %f; data cost = %f; reg cost = %f. \n',d+e,rms,d,e);
    
    theta=W.*W;
    
    tmp=ft'*(ftt'*v);
    tmp=tmp./(lambda+theta);
    rho=theta.*tmp;
    
end

res=ftt*rho+mres;

figure;
imshow(abs(res(:,:,1,1)),[]);

disp('');


function cov=get_cov(res)

res=abs(res);
nt=size(res,4);

ma=mean(res,4);
sd=std(res,[],4);
sd2=sd.^2;
cov=zeros(nt,nt);
m=abs(ma)>max(abs(ma(:)))*0.1;
for i=1:nt
    for j=i:nt
         tmp=(res(:,:,:,i)-ma).*(res(:,:,:,j)-ma)./sd2;
         cov(i,j)=mean(tmp(m));
         cov(j,i)=cov(i,j);             
    end
end
        
        
        
        