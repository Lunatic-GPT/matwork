function ts=peak_remove(ts,f,sw)

% d=peak_remove(ts,f,sw)

% f is the frequency range of the peak. 0 at center
% sw is the total spectral width.
s=fft(ts);

s=fftshift(s);

N = length(ts);
if mod(N,2)==0
   ind =-N/2+1:N/2;
else
    ind = -(N-1)/2:(N-1)/2;
end


fa=ind*sw/N;

ind2=fa>f(1) & fa<f(2);

in1=max(find(ind2))+1;
in2=min(find(ind2))-1;
s(ind2)=(s(in1)+s(in2))/2;
s=fftshift(s);

ts=ifft(s);


    


    




