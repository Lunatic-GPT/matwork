function recon_fl_fq_retro_async(dname)

load(dname);
nch=length(unique(Channel));
physio=icePara(1:nch:end,4);

minInterval=0.4;
TR=0.03;
nphase=10;
trig_delay=0;
thr=0;
plot_results=1;

[ph,trig,hr]=physioTrigger(physio,minInterval,TR,nphase,trig_delay,thr,plot_results);

Line1c=Line(1:nch:end);
nLine=length(unique(Line1c));
m=zeros(nphase,nLine);

for i=1:length(ph)
   m(ph(i), Line1c(i)+1)=m(ph(i), Line1c(i)+1)+1;   
end
figure;imshow(m==0,[]);

figure;imshow(m,[]);
disp('');
