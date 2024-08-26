function fit_mTE(prefix,tind,nte_exclude)
% t2(prefix,tind[,nte_exclude])

if ~exist('nte_exclude','var')
    nte_exclude=[];
end

fname = [prefix,'.fid'];
a=read_fid(fullfile(fname,'fid'));
a=squeeze(a);
a=a-repmat(mean(a(end-100:end,:),1),[size(a,1),1]);
te = readPar(fname,'te');
aa = abs(a);
b= mean(aa(tind,:),1);
figure;
plot(real(a(:,1)),'r');
hold on;
plot(imag(a(:,1)),'b');
plot(abs(a(:,1)),'k');

plot([min(tind)+1,min(tind)+1],ylim,'k');
plot([max(tind)+1,max(tind)+1],ylim,'k');

figure;
plot(te,log10(b),'o');
[te,te_ind] = sort(te);
b = b(te_ind);
beta=nlinfit(te(1:end-nte_exclude+1),b(1:end-nte_exclude+1),@exp_decay,[max(b),0.1]);
hold on;
fprintf('T2 = %4.3f s\n',beta(2));
x=linspace(0,max(te));
plot(x,log10(exp_decay(beta,x)),'k');




