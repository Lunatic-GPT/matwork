% refer to the circuit diagram on Page 70.
%c1 match, c3 tuning.
c2=1e-12;
c4=2e-12;
%c2=1e-14;
%c4=1;
%c3=2.5e-12;
%c1=linspace(0,15e-12,100);
c3=linspace(0,15e-12);
%[c1,c3]=meshgrid(x,y);
w = 400.0e6;
L = 9.6e-7;%1/w/w/(c3+c2*c4/(c2+c4));
r=2;
z_match=zeros(size(c1));
c1_match=zeros(size(c1));
z2=zeros(size(c3));
for i=1:length(c3)
      
        z_tmp=1/(1i*w*c4)+z_par(1/(1i*w*c3(i)),1i*w*L+r);
        z_tmp2=z_par(1/(1i*w*c2),z_tmp);
        
        c1=1/(imag(z_tmp2)*w);
        z2(i)=real(1/(1i*w*c1)+z_tmp2);
        c1_match(i)=c1;           
end

c1_match(z2>40 & z2<60)

%figure;plot(c3,z2);

%figure;plot(c3,c1_match);
