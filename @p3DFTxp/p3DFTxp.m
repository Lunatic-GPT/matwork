function  res = p3DFTxp(mask,phase,mode)

%res = p2DFT(mask[,phase,mode])
%
%
%	Implementation of partial Fourier operator.
%	
if nargin <2
	phase = 1;
end
if nargin <3
	mode = 2; % 0 - positive, 1- real, 3-cmplx
end


res.adjoint = 0;
res.mask = mask;
res.imSize = size(mask);
res.dataSize = size(mask);
res.ph = phase;
res.mode = mode;
res = class(res,'p3DFTxp');

