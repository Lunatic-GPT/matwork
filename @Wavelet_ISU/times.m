function res = times(a,b)



if isa(a,'Wavelet_ISU') == 0
    error('In  A.*B only A can be Wavelet operator');
end

if a.adjoint
    res = a.wav'*b*a.wav;
else
    res = a.wav*(b)*a.wav';
end


