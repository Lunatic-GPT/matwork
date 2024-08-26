% sinc pulse simulation;

gamma = 2*pi*4.257; %rad.kHz/Gauss
t = 0.003; % seconds
N = 250;
rf = 0.0782916*msinc(N,1)*3;  % Gauss
gz = 1*ones(size(rf)); % Gauss/cm
T1=1.01;
T2 = 0.2;
z = -1:0.01:1;  % cm
[mx,my,mz] = bloch(rf,gz,t/N,T1,T2,0,z,0);

figure;plot(z,[mx(:),my(:),mz(:)]);
legend('mx','my','mz');