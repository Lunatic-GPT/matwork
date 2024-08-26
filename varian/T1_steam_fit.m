function T1_steam_fit(fname,tind)
% T1_steam_fit(fname,tind)
% tind: the time points used to calculate the recovered magnetization.


if nargin ==0
    help T1_steam_fit
    return;
end
a=read_fid(fullfile(fname,'fid'));
a=squeeze(a);
a=a-repmat(mean(a(end-100:end,:),1),[size(a,1),1]);
ti = readPar(fullfile(fname,'procpar'),'ti');

[tmp,ind]=max(ti);
ph=angle(a(tind(1),ind));
a=a*exp(-1i*ph);
aa=real(a);


b= mean(aa(tind,:),1);


figure;
plot(ti,b,'o');
beta=nlinfit(ti,b,@T1_rcvr,[max(b),1,2]);
fprintf('Sat level = %4.2f%%, T1 = %4.3f\n',beta(2)*100,beta(3));
hold on;
x=linspace(0,max(ti),100);
plot(x,T1_rcvr(beta,x),'r-');

[tmp,i_maxti]=max(ti);
b2=1-b/b(i_maxti);
figure;
plot(ti,log10(b2),'o');
hold on;
yfit = log10(beta(2)*exp(-ti/beta(3)));
plot(ti,yfit,'k-');


