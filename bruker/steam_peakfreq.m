function steam_peakfreq(prefix)
a=readb_fid(prefix);
dly=readbPar([prefix,'/method'],'PVM_DigShift');
a=a(:);
a=a(dly+1:end);
sw=readbPar([prefix,'/method'],'PVM_DigSw');

a=a*exp(-1i*148/180*pi);

fa=fft1c(a,1);

[tmp,imax]=max(abs(fa));

ph=plot_sp(a,0,sw,false);


bf=readbPar([prefix,'/acqp'],'BF1');

f=bf+(imax-length(a)/2-1)*sw/length(a)*1.0e-6;

fprintf('%9.6f\n',f);





