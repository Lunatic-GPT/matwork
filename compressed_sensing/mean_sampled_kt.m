function res=mean_sampled_kt(z,m)

nd=ndims(m);
z=z.*m;
mm=sum(m,nd);

mz=sum(z,nd);
mm(mm==0)=1;

res=mz./mm;



