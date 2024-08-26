function [mx,my,mz]=rf_simulate_siemens(fname,pw,B1max,freq,M0)
% rf_simulate_varian(fname,pw,B1max,freq,M0)
% fname: name of the rf pattern
% pw : total pulse width in seconds.
% B1max: peak RF field strength in Hz. (Larmor frequency of proton at that field) 
% freq: frequencies to simulate. in Hz.
% M0: a 1*3 matrix specifying the initial magnetization. default:[0,0,1]



T1 = 2.01;
T2 = 0.2;


B1max = B1max/4258;  %convert to Gauss.

rf=read_rf_siemens(fname);
rf=rf*B1max;

t=pw*ones(size(rf,1),1)/size(rf,1);

if ~exist('M0','var')
   [mx,my,mz]=bloch(rf,zeros(size(rf)),t,T1,T2,freq(:),0,0);
else
    M0 = repmat(M0(:),[1,length(freq)]);
   [mx,my,mz]=bloch(rf,zeros(size(rf)),t,T1,T2,freq(:),0,0,M0(1,:),M0(2,:),M0(3,:));
end    
figure;subplot(2,1,1);
plot(freq,mz(:),'k','LineWidth',2);
hold on;plot(freq,sqrt( my(:).^2+mx(:).^2),'r','LineWidth',2);
legend('mz','mxy');
xlabel('Frequency (Hz)');

ph = atan2(my,mx)*180/pi;
subplot(2,1,2);plot(freq,ph,'g','LineWidth',1);
xlabel('Frequency (Hz)');
ylabel('phase');
ylim([-180,180]);


