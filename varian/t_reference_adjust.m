function t_reference_adjust(h)

h_toff=findobj(h,'Tag','toff');
h_toff_edit=findobj(h,'Tag','toff_edit');

h_sigma_edit=findobj(h,'Tag','sigma_edit');
sigma=get(h_sigma_edit,'String');
sigma=str2double(sigma);

if h_toff==gcbo
    toff = get(gcbo,'Value');
else
    toff=get(gcbo,'String');
    toff=str2double(toff);
end

set(h_toff,'Value',toff);
set(h_toff_edit,'Value',toff,'String',sprintf('%d',toff));

hr = findobj(h,'Tag','fid_real');
hi = findobj(h,'Tag','fid_imaginary');
r = get(hr,'YData');
im = get(hi,'YData');
x = get(hr,'XData');

h_ref = findobj(h,'Tag','fid_reference');
set(h_ref,'XData',[toff+1,toff+1]*(x(2)-x(1)));


d = r(toff+1:end)+1i*im(toff+1:end);
d(end+1:end+toff) = 0;
%r=circshift(r,[1,-toff]);
%im=circshift(im,[1,-toff]);
t=1:length(d);
d2=d.*exp(-((t-1)/sigma).^2);
f=fft(d2);
f=fftshift(f);

% recalculate spectrum
hr = findobj(h,'Tag','sreal');
hi = findobj(h,'Tag','simaginary');
set(hr,'YData',real(f));
set(hi,'YData',imag(f));


