function [theta,phi]=unitVec2thetaPhi(x)
%  [theta,phi]=unitVec2thetaPhi(x)
%theta and phi in degree;
% x is a 3*n matrix

x=x./repmat(sos(x,1),[size(x,1),1]);

r=sos(x(1:2,:),1);

theta=real(asin(r));

theta(x(3,:)<0)=pi-theta(x(3,:)<0);

phi=atan2(x(2,:),x(1,:));

theta=theta*180/pi;

phi(~any(abs(x(1:2,:))>1e-6,1))=0;

phi=phi*180/pi;



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
