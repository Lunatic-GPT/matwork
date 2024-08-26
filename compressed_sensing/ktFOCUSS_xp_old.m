function ktFOCUSS_xp(v,m,lambda,p,rho)

if ~exist('p','var')
    p=1/2;
end


sz=size(m);
ft=p2DFTxp(m,sz(1:2),1,2);
ftt=FTt;


if ~exist('rho','var')
   v=v.*m;
   mv=sum(v,4)./sum(m,4);
   rho = ft'*(ftt'*mv); 
end

for i=1:10
   
    W=abs(rho).^p;    
    
    d=(v-ft*(ftt*rho)).^2;
    d=sum(d(:));
    e=lambda*rho./W;
    e=sum(e(:));
  
    fprintf('Cost = %f\n',d+e);
    
    theta=W*W';
    
    tmp=ft'*(ftt'*v);
    tmp=tmp./(lambda+theta);
    rho=theta.*tmp;
    
end
