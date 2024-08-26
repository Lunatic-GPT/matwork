function mask=cstable_1dpower(npe,reduction,etl,alpha,fname)
%mask=cstable_1dpower(npe,reduction,etl,alpha,fname)
% generate random k-space lines based on a power distribution 1/x^alpha.

fid=fopen(fname,'w');
fprintf(fid,'t1 = \n');


petable=[];


c=fzero(@(x) area_power(x,alpha,npe)-1/reduction,1);

while length(petable)<npe/reduction
    
    x=linspace(-npe/2,npe/2-1,npe);
      
    y=c./abs(x).^alpha;
    y(npe/2+1)=1;
    
    
    mask=rand(1,npe)<y;
  %  mask=rand(1,npe)<1/reduction;
    tmp=find(mask>0);
    disp(length(tmp));
    tmp2=setdiff(tmp,petable);
    petable=[petable,tmp2];
    
end

petable=petable(1:npe/reduction)-npe/2;


[tmp,ind]=sort(abs(petable));

nseg=npe/reduction/etl;

y=[];
for i=1:nseg
   for j=1:etl             
       ii=ind((j-1)*nseg+i);
      fprintf(fid,'  %d\n',petable(ii));
      y(end+1)=petable(ii);
   end
end

fclose(fid);

mask=zeros(1,npe);
mask(petable+npe/2)=1;
mask=repmat(mask,[npe,1]);
figure;imagesc(mask);




function res=area_power(c,alpha,n)

a=linspace(-n/2,n/2-1,n);

b=c./abs(a).^alpha;

b(n/2+1)=1;
b(b>1)=1;

res=sum(b)/n;




