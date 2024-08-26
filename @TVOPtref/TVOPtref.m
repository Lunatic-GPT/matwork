function  res = TVOPtref(ref)

%res = TVOP()
%
% Implements a spatial finite-differencing operator.
%
% (c) Michael Lustig 2007

res.adjoint = 0;
res.ref=ref;
res = class(res,'TVOPtref');

