function [mask,b]=cstable_2dpower(npe,reduction,alpha,etl,fname)

% generate random k-space lines based on a gaussian distribution.
%mask=cstable_2dpower(npe,reduction,alpha,etl,fname)


petable=[];


c=fzero(@(x) area_power(x,alpha,npe)-1/reduction,39/reduction);

while length(petable)<npe*npe/reduction
    
    x=linspace(-npe/2,npe/2-1,npe);
      
    y=c./abs(x).^alpha;
    y(npe/2+1)=1;
    
     a=linspace(-npe/2,npe/2-1,npe);
     [x,y]=meshgrid(a);

     b=c./sqrt(x.^2+y.^2)^alpha;

     b(npe/2+1,npe/2+1)=1;
     b(b>1)=1;
    mask=rand(npe,npe)<b;
  %  mask=rand(1,npe)<1/reduction;
    tmp=find(mask>0);
    tmp2=setdiff(tmp,petable);
    petable=[petable;tmp2];
    
end


[tmp,ind]=sort(abs(petable-npe*npe/2));

nseg=npe*npe/reduction/etl;

if exist('fname','var')
fid=fopen(fname,'w');
for i=1:nseg
   for j=1:etl             
       ii=ind((j-1)*nseg+i);
      fprintf(fid,'  %d\n',petable(ii));
   end
end

fclose(fid);

end

m=zeros(npe,npe);
petable=petable(1:npe*npe/reduction);
m(petable)=1;


%figure;imagesc(m);
%xlim([0.5,npe+0.5]);
%ylim([0.5,npe+0.5]);






function res=area_power(c,alpha,n)

a=linspace(-n/2,n/2-1,n);
[x,y]=meshgrid(a);


b=c./sqrt(x.^2+y.^2)^alpha;
     
b(n/2+1,n/2+1)=1;
b(b>1)=1;

res=sum(b(:))/n/n;




