function sat_rcvr_fit_fid(fid_prefix,tind)
% sat_rcvr_fit(fid_prefix,tind)
a=read_fid([fid_prefix,'.fid/fid']);
ti=readPar(fid_prefix,'ti');
bl=mean(a(end-99:end,1,:),1);
a=a-repmat(bl,[size(a,1),1,1]);
aa=abs(a);
area=sum(aa(tind,1,:),1);
plot(ti,squeeze(area),'o');
b=nlinfit(ti(:),area(:),@T1_rcvr,[3000000,1,2]);
hold on;
fprintf('sat = %4.3f, T1 = %4.3s s\n',b(2),b(3));

x=linspace(0,max(ti),100);
plot(x,T1_rcvr(b,x),'r-');