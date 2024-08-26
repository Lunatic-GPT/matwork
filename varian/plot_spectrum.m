function plot_spectrum(fid,index,toff,phase)
% plot_spectrum(fid,index,toff,phase)

if ~exist('index','var')
    index=1;
end
a=read_fid(fullfile(fid,'fid'));
 sw=readPar(fid,'sw');
 a=squeeze(a);
 a=a(:,index);

a=a-mean(a(end-99:end));

if ~exist('toff','var')
    
    [tmp,toff] = max(abs(a));
    
    toff = toff -1;
    phase = phase_adjust_auto(a(toff+1:end));
end

 a= sp_phasecorr(a,phase);
 dt=1000/sw;
t = (1:size(a,1))*dt;
sigma = max(t)/5;

 figure;
 set(gcf,'ToolBar','figure');
 set(gcf,'Position',[393   185   560   600]);
 set(gcf,'UserData',phase);
 
h_phase_edit = uicontrol(gcf,'Style','edit','Value',phase,'Tag','phase_edit','String',num2str(phase));
set(h_phase_edit,'Callback',sprintf('phase_adjust(%d)',gcf));
h_phase=uicontrol(gcf,'Style','slider','Max',180,'Min',-180,'Tag','phase','Value',phase);
set(h_phase,'Position',[80,20,150,20],'SliderStep',[1/360,1/36]);
set(h_phase,'Callback',sprintf('phase_adjust(%d)',gcf));

h_toff_edit = uicontrol(gcf,'Style','edit','Value',toff,'Tag','toff_edit','String',num2str(toff));
set(h_toff_edit,'Position',[290,20,40,20]);
set(h_toff_edit,'Callback',sprintf('t_reference_adjust(%d)',gcf));
h_toff=uicontrol(gcf,'Style','slider','Max',size(a,1),'Min',0,'Tag','toff','Value',toff);
set(h_toff,'Position',[340,20,150,20],'SliderStep',[1/size(a,1),10/size(a,1)]);
set(h_toff,'Callback',sprintf('t_reference_adjust(%d)',gcf));

h_sigma_edit = uicontrol(gcf,'Style','edit','Value',sigma,'Tag','sigma_edit','String',num2str(sigma));
set(h_sigma_edit,'Position',[80,50,40,20]);
set(h_sigma_edit,'Callback',sprintf('sigma_adjust(%d)',gcf));
h_sigma=uicontrol(gcf,'Style','slider','Max',size(a,1)*dt,'Min',0,'Tag','sigma','Value',sigma);
set(h_sigma,'Position',[130,50,200,20],'SliderStep',[1/size(a,1),10/size(a,1)]);
set(h_sigma,'Callback',sprintf('sigma_adjust(%d)',gcf));

subplot(2,1,1);

plot(t,real(a),'r','Tag','fid_real');
xlabel('Time (ms)');
hold on;
plot(t,imag(a),'b','Tag','fid_imaginary');
%plot(t,abs(a),'k','Tag','fid_abs');
yl=ylim;
plot([t(1),t(1)]*(toff+1),yl,'k','Tag','fid_reference');

set(gca,'Position',[0.1300    0.6499    0.7750    0.3]);


xlim([0,max(t)]);
legend('real','imag','abs');

a2 = a(toff+1:end);
a2(end+1:end+toff)=0;

nt2=length(a2);
t2=(1:nt2)*dt;

yl=ylim;
plot(t2,yl(2)/2*exp(-((t2-dt)/sigma).^2),'k-','Tag','Gauss');
t2=reshape(t2,size(a2));
a2=a2.*exp(-((t2-dt)/sigma).^2);
%a2 = circshift(a,-toff);
f=fft(a2); %fft(a2) and fft(a2') give reversed results.

f=fftshift(f);

subplot(2,1,2);

xf = (-nt2/2:nt2/2-1)/nt2*sw;
xf = -xf/400+4.75;  %assuming water is at resonance.

plot(xf,imag(f),'w','Tag','simaginary');
hold on;
plot(xf,real(f),'r','Tag','sreal');

set(gca,'XDir','reverse');
set(gca,'Position',[0.1300    0.2499    0.7750    0.3]);
xlim([min(xf),max(xf)]);
xlabel('Chemical Shift (ppm)');



   