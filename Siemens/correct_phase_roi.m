function correct_phase_roi(dname,roi_name)

m=load(roi_name);

m=clusterize2(m.roi);

a=ri(dname);

a2=mean(a,4);

x=[];
y=[];

for i=1:max(m(:))
    
pos=ind2subb(size(m),find(m==i));
x(i,:)=mean(pos,1);

y(i)=mean(a2(m==i));

end

x=[x,ones(size(x,1),1)];

b=x\y';

%%

figure;plot(y,x*b,'o');


for k=1:size(a,4)
  for i=1:size(a,1)
    for j=1:size(a,2)     
        a(i,j,1,k)=a(i,j,1,k)-[i,j,1]*b;
    end
  end
end
dname=strtok(dname,'.');

writeanalyze(a,[dname,'_detrend']);

