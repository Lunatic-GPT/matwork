function spa_mrui(fname,frange,fshift,fid_prefix)
% spa_mrui(fname,frange,fshift,fid_prefix)
% plot an spectrum array in mrui format.
% fname: name of the mrui spectrum array
% frange: the frequency range to plot
% fshift: right shift of consecutive spectrum
% fid_prefix: used to determine the carrier frquency.
%Lac 1.5
% NAA: 2.1
% Glu: 3.78.
% Cr: 3.06, 3.95

satof = readPar(fid_prefix,'satof');
tof = readPar(fid_prefix,'tof');
f0=(tof-satof)/400+4.7;

[h,a]=read_mrui(fname);

sw=10;  % 10 ppm
nt=size(a,2);

a=a-repmat(mean(a(1600:end,:),1),[2043,1]);

t=(1:2043)'/4000;
a=a.*repmat(exp(-t.^2/2/0.2^2),[1,nt]);
fa=fft(a,[],1);
fa=fftshift(fa,1);
%%
xf=sw*((1:2043)-1022)/2043+f0;

ind=find(frange(1)<xf&frange(2)>xf);



figure;

for i=1:size(a,2)
    
plot(xf(ind)-(i-1)*fshift,real(fa(ind,i)),'r-');
hold on;
end
set(gca,'XDir', 'reverse','FontSize',16,'XTick',0:5);
xlim([frange(1)-fshift*(size(a,2)-1),frange(2)]);
%ylim([-0.1,0.5]);

