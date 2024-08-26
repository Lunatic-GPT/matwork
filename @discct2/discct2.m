function res = discct2(mn)
% res = Wavelet(Filtertype, filterSize, wavScale)
%
% implements a wavelet operator
%
% (c) Michael Lustig 2007

res.adjoint = 0;
res.m=mn(1);
res.n=mn(2);
res = class(res,'discct2');
