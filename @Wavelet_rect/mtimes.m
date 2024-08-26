function res = mtimes(a,b)

if isa(a,'Wavelet_rect') == 0
    error('In  A.*B only A can be Wavelet operator');
end

if a.adjoint
    res = IWT2_PO_rect(real(b),a.wavScale,a.qmf) + 1i*IWT2_PO_rect(imag(b),a.wavScale,a.qmf);
else
    res = FWT2_PO_rect(real(b),a.wavScale,a.qmf) + 1i* FWT2_PO_rect(imag(b),a.wavScale,a.qmf);
end


