function y=hypergeometric2F1(a,b,c,z)

if abs(z) >=1
    error('not converge for abs(z)>=1');
end

y=1;

if z==0
    return;
end


k = 0;
tmp =1;
while 1
    
   tmp = tmp*(a+k)*(b+k)*z/(c+k)/(k+1);
   tmp_sum = tmp;
    k = k+1;
   tmp = tmp*(a+k)*(b+k)*z/(c+k)/(k+1);
   tmp_sum = tmp_sum+tmp;
   k = k+1;
   if abs(tmp_sum)< 1/100000
       break;
   end
   y = y+tmp_sum;
  
end
       



