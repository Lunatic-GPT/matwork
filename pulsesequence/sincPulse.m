function y=sincPulse(pw,nc,N,B1max)
% y=sincPulse(pw,nc,N,B1max)
% pw: pulse width.
% nc: number of cycles on one side (positive or negative).
% N total number of time points.
% B1max: unit Hz


x = linspace(-pw/2,pw/2,N);

df = nc*4/pw;

fprintf('band width = %4.3f Hz\n',df);
y = sinc(x*df)*B1max/4258;


figure;
plot(x,y);
xlabel('Time (s)');
ylabel('RF field (Gauss)');

T1 = 1.01;
T2 = 0.2;
freq = linspace(-3*df,3*df,100);

[mx,my,mz]=bloch(y,zeros(size(y)),pw/N,T1,T2,freq,0,0,0,0,1);

figure;plot(freq,[mx(:) my(:) mz(:)]);
legend('mx','my','mz');
xlabel('Frequency (Hz)');