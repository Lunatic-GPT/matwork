sz = 128;
k=0;
ProUndersampling=16;
  angle = 2 * pi;
  n_theta = round(pi*sz/(2*sqrt(ProUndersampling)));
x=[];
y=[];
z=[];
nphi=[];
thetaa=[];
for j=1:n_theta-1
  
    n_phi = round(pi*sz*sin(pi*j/n_theta))/sqrt(ProUndersampling);

    theta=[];
    for i=0:n_phi-1
    
      x(end+1)=sin(pi*j/n_theta)*cos(angle*i/n_phi);
      y(end+1)=sin(pi*j/n_theta)*sin(angle*i/n_phi);
      z(end+1)=cos(pi*j/n_theta);
      sinx=z(end)./sqrt(x(end).^2+y(end).^2+z(end).^2);
      theta(end+1)=asin(sinx);
      

    end
    thetaa(end+1)=theta(1);
 %%  figure; plot(theta);
 nphi(end+1)=n_phi;    
end

phi=atan2(x,y);
figure;hist(phi);
sinx=z./sqrt(x.^2+y.^2+z.^2);
theta=asin(sinx);
figure;hist(theta,-pi/2:0.0628:pi/2);
