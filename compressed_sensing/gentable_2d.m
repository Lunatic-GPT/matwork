function gentable_2d(npe,reduction,etl,fname)

% generate random k-space lines based on a gaussian distribution.

fid=fopen(fname,'w');
fprintf(fid,'t1 = \n');


petable=[];


sig=fsolve(@(x) area_gauss(x,npe)-1/reduction,0.5);

while length(petable)<npe/reduction
    
    x=(-npe/2:npe/2-1)/npe*2;
  
  mask=rand(1,npe)<exp(-x.^2/2/sig/sig);
  %  mask=rand(1,npe)<1/reduction;
    tmp=find(mask>0);
    disp(length(tmp));
    tmp2=setdiff(tmp,petable+npe/2);
    
    petable=[petable,tmp2-npe/2];
    
end

petable=petable(1:npe/reduction);
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

mask=zeros(1,npe);
mask(petable+npe/2)=1;
mask=repmat(mask,[npe,1]);
figure;imagesc(mask);

fclose(fid);



function frac=area_gauss(sig,n)

a=linspace(-n/2,n/2-1,n);
a=a/n*2;

res=sum(exp(-a.^2/sig/sig/2));
frac=res/n;
