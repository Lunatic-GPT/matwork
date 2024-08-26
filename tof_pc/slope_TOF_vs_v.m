function k=slope_TOF_vs_v(thk,TR,fa,T1)
v=linspace(0.01,0.3,20);
nv=length(v);

v=[v,linspace(0.3,1.5,20)];

for i=1:length(v)
  sart(i)= TOF_signalIntensity(v(i),thk,TR,fa,T1);
end

xmat=[ones(nv,1),v(1:nv)'];

b=xmat\sart(1:nv)';

figure;plot(v,sart,'o-');
hold on;


xmat=[ones(length(v),1),v(:)];
plot(v,xmat*b,'r-');
k=b(2)/b(1);