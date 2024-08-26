function res = fft2o()
% res = Wavelet(Filtertype, filterSize, wavScale)
%
% implements a wavelet operator
%
% (c) Michael Lustig 2007

res.adjoint = 0;
res = class(res,'fft2o');
