function d2=meanc(d,c,vargin)


for i=1:max(c)
   
d2(i)=mean(d(c==i),vargin);    
end

