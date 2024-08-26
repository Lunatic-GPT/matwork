function pm=Daubechies4_1D(m)
       
if size(m,1)==1
    col=false;
    m=m(:);
else
    col=true;
end


       N = length(m);
       s1 = m(1:2:N-1) + sqrt(3)*m(2:2:N);
       d1 = m(2:2:N) - sqrt(3)/4*s1 - (sqrt(3)-2)/4*[s1(N/2); s1(1:N/2-1)];
       s2 = s1 - [d1(2:N/2); d1(1)];
       s = (sqrt(3)-1)/sqrt(2) * s2;
       d = (sqrt(3)+1)/sqrt(2) * d1;
       pm=[s;d];
       
       if ~col
           pm=pm';
       end

