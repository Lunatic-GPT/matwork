function correct_intensity_inhomogeneity(fname, m,order)

if isa(m,'char')
m=ri(m);
end

m=clusterize2(m);

a=ri(fname);

x=[];
y=[];

for i=1:max(m(:))
    
pos=ind2subb(size(m),find(m==i));
x(i,:)=mean(pos,1);

y(i)=mean(a(m==i));

end

xv=x;


x=[xv,ones(size(xv,1),1)];


if order>1    
    x=[x,xv(:,1).^2,xv(:,2).^2,xv(:,1).*xv(:,2)];
end

    
mi=mean_roi(a,m>0);

b=x\y';

a2=0*a;
for i=1:size(m,1)
    for j=1:size(m,2)
   
        
        if order ==1
         a2(i,j)=a(i,j)-[i,j,1]*b+mi;
        else
         a2(i,j)=a(i,j)-[i,j,1,i*i,j*j,i*j]*b+mi;
        end
    end
    
end

fname=strtok(fname,'.');
save([fname,'_corr'],'a2');





