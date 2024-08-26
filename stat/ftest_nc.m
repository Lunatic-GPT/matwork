function p = ftest_nc(v1,v2,x,ld)
% p = ftest_nc(v1,v2,x,ld)

p = zeros(size(x));

for i=1:length(x)
  p(i) = 1-quad(@(z) fdist_nc(z,v1,v2,ld) ,0,x(i));
end
    
    
    

