function b=fid_undersample(a,factor,np)
% a: the origina data
% undersampling factor: should be greater than 1.
% np: number of data points in the output.
b=zeros(1,np);
b(1)=a(1);

for i=1:np-1
   
  ind=round(factor*i);  
  if abs(round(ind)-factor*i)<1e-5
     b(i+1)=a(ind+1);
  else
     l1=ceil(factor*i);
     l2=floor(factor*i);
     
     b(i+1)=a(l2+1)*(l1-factor*i)+a(l1+1)*(factor*i-l2);
     
  end
    
end