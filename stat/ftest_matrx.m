function [p,f]=ftest_matrx(data)
% [p,f]=ftest_matrx(data)
% data is a 2-dimensional matrix where each column is a separate treatment.

sz1=size(data,1);
sz2=size(data,2);
s_wit = zeros(1,sz2);
dmean = mean(data,1);
for i = 1:sz2
    
    s_wit(i) = std(data(:,i));
end

   f = sz1*std(dmean)^2/mean(s_wit.^2);
   if f==0 || isnan(f) 
       p = 1;
       f=0;
       return;
   end
   dof1 = sz2-1;
   dof2 = sz2*(sz1-1);
   p = FTest(dof1,dof2,f);
   