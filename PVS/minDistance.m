function res=minDistance(x0,v,x1)
% find the minimum distrance between v and a line that is in the direction of v and passes through x0.

% the line is x=x0+v*t;
v=v/sos(v);
tmin=sum(v.*(x1-x0));
res=sos(x0+v*tmin-x1);


function res = sos(x ,dim, pnorm)
% res = sos(x [,dim, pnorm])
%
% function computes the square root of sum of squares along dimension dim.
% If dim is not specified, it computes it along the last dimension.
%
% (c) Michael Lustig 2009

if nargin < 2
    dim = size(size(x),2);
end

if nargin < 3
    pnorm = 2;
end
res = (sum(abs(x.^pnorm),dim)).^(1/pnorm);
