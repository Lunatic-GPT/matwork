function sigma_adjust(h)

h_toff=findobj(h,'Tag','toff');
toff = get(h_toff,'Value');
    
h_sigma_edit=findobj(h,'Tag','sigma_edit');
h_sigma=findobj(h,'Tag','sigma');

if h_sigma==gcbo
    sigma = get(gcbo,'Value');
else  
  sigma=get(h_sigma_edit,'String');
  sigma=str2double(sigma);
end

set(h_sigma,'Value',sigma);
set(h_sigma_edit,'Value',sigma,'String',sprintf('%f',sigma));

hr = findobj(h,'Tag','fid_real');
hi = findobj(h,'Tag','fid_imaginary');
r = get(hr,'YData');
im = get(hi,'YData');
x = get(hi,'XData');
dt=x(2)-x(1);

hs = findobj(h,'Tag','Gauss');
t2=(1:length(r));
ax=get(hr,'Parent');
axes(ax);
yl=ylim;
y = yl(2)/2*exp(-((t2-toff-1)*dt/sigma).^2);
set(hs,'YData',y);
% recalculate spectrum

d = r(toff+1:end)+1i*im(toff+1:end);
d(end+1:end+toff) = 0;
%r=circshift(r,[1,-toff]);
%im=circshift(im,[1,-toff]);
t=1:length(d);
%d2=d;
d2=d.*exp(-((t-1)*dt/sigma).^2);
f=fft(d2);
f=fftshift(f);
hr = findobj(h,'Tag','sreal');
hi = findobj(h,'Tag','simaginary');
set(hr,'YData',real(f));
set(hi,'YData',imag(f));


