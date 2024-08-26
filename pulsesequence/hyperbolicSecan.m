function [m,ph]=hyperbolicSecan(tdf,df,cutoff,B1max,npts)
% hyperbolicSecan(tdf,t,cutoff,B1max,npts)
% tdf: Time * bandwidth (typical value: 10)
% df : bandwidth in Hz
% cutoff: cutoff x value of sech(x). (typical value: 6)
% B1max : maximum field in Hz
% npts: number of points

t = tdf/df;
beta = 2*cutoff/t;

mu = df*pi/beta;
A = B1max/4258;

A2 = sqrt(mu)*beta/2/pi;

fprintf('B1max = %f (should be >> %f)\n',B1max,A2);

ta = (-(npts-1)/2:(npts-1)/2)*t/npts;

m = A*sech(beta*ta);

ph = mu*log(sech(beta*ta));

figure;plot(ta,m*4258);
title('RF strength');
xlabel('time (s)');
ylabel('Hz (strength)');

figure;plot(ta,ph);
title('RF Phase');
xlabel('time (s)');
ylabel('phase (rad)');

rf = m.*cos(ph)+1i*m.*sin(ph);

T1 = 1.01;
T2 = 0.2;
freq = linspace(-3*df,3*df,100);
[mx,my,mz]=bloch(rf,zeros(size(rf)),t/npts,T1,T2,freq,0,0,0,0,1);

figure;plot(freq,[mx(:) my(:) mz(:)]);
legend('mx','my','mz');
xlabel('Frequency (Hz)');



