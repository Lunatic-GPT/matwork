function PwrCal(fname)
a=readb_fid(fname);

pref=readbPar(fullfile(fname,'method'),'PVM_RefAttCh1');

pwr=readbPar(fullfile(fname,'method'),'pwrlist');

pwr=pref-pwr;

pwr=pwr';
x=db2factor(pwr);

%y=squeeze(abs(a(31,1,:)));
%y=y';
fa=fft1c(a,1);

%y=abs(fa(31,:));

y=mean(squeeze(abs(fa)),1);

figure;plot(x,y,'ko');


options=optimset('MaxFunEvals',40);
beta=lsqcurvefit(@sinfita,[max(y),pi/2],x,y,[],[],options);
disp(-factor2db(beta(2)/pi*2)+pref);
hold on;
plot(x,sinfita(beta,x),'r');

