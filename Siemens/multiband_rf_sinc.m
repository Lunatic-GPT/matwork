function out=multiband_rf_sinc(n,df,ph)
% n: number of bands
% df: frequency seperations (in Hz) between bands for 1 ms pulse
% The single band sinc pulse has a frequence width of 6.2 kHz for a 1 ms pulse.
% ph: phase for each band; 1*n; unit: degree

 a=read_rf_bruker('e:/Dropbox/ForBruker/wave/sinc.exc',1);
 
 ib=-(n-1)/2:(n-1)/2;
     
 b=zeros(size(a,1),n);
 
     for i=1:n
         
         lf=df/1000*ib(i)*2*pi*(0:length(a)-1)/length(a);
         b(:,i)=a.*exp(1i*lf')*exp(1i*ph(i)/180*pi);
     end
     
out=sum(b,2);

