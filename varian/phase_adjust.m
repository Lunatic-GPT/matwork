function phase_adjust(h)

h_phase=findobj(h,'Tag','phase');
h_phase_edit=findobj(h,'Tag','phase_edit');

if h_phase==gcbo
    phase = get(gcbo,'Value');
else
    phase=get(gcbo,'String');
    phase=str2double(phase);
end

set(h_phase,'Value',phase);
set(h_phase_edit,'Value',phase,'String',sprintf('%4.1f',phase));

% adjust fid
hr = findobj(h,'Tag','fid_real');
hi = findobj(h,'Tag','fid_imaginary');
r = get(hr,'YData');
im = get(hi,'YData');

phase_current = get(h,'UserData');

dph = phase-phase_current;

d=(r+1i*im);
a = sp_phasecorr(d,dph,0);
set(hr,'YData',real(a));
set(hi,'YData',imag(a));
set(h,'UserData',phase);

% adjust spectrum
hr = findobj(h,'Tag','sreal');
hi = findobj(h,'Tag','simaginary');
r = get(hr,'YData');
im = get(hi,'YData');
d=(r+1i*im);
a = sp_phasecorr(d,dph,0);
set(hr,'YData',real(a));
set(hi,'YData',imag(a));



