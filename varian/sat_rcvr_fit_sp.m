function sat_rcvr_fit_sp(fname,toff,f_range,phase)
% T1_sat_rcvr_fit(fname,t_off,f_range,phase)
% t_off: the first time points to perform fourier transformation.


if nargin ==0
    help T1_steam_fit
    return;
end
a=read_fid(fullfile(fname,'fid'));
a=squeeze(a);
a=a-repmat(mean(a(end-100:end,:),1),[size(a,1),1]);

a= a*exp(1i*phase*pi/180);
 
ti = readPar(fname,'ti');
 sw=readPar(fname,'sw');
a2 = a(toff+1:end,:);
a2(end+1:end+toff,:)=0;
f=fft(a2,[],1);
f=fftshift(f,1);


xf = (-size(a,1)/2:size(a,1)/2-1)/size(a,1)*sw;
xf = xf/400+4.8;  %assuming water is at resonance.

rf=real(f);
f_ind = find(xf>=f_range(1)&xf<f_range(2));
b= mean(rf(f_ind,:),1);


figure;
plot(ti,b,'o');
beta=nlinfit(ti,b,@T1_rcvr,[max(b),1,median(ti)]);
fprintf('Sat level = %4.2f%%, T1 = %4.3f s\n',beta(2)*100,beta(3));
hold on;
x=linspace(0,max(ti),100);
plot(x,T1_rcvr(beta,x),'r-');

[tmp,i_maxti]=max(ti);
%b2=1-b/beta(1);
b2=1-b/b(i_maxti);

figure;
plot(ti,log10(b2),'o');
hold on;
yfit = log10(beta(2)*exp(-ti/beta(3)));
plot(ti,yfit,'k-');


