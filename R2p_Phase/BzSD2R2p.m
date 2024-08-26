function R2p=BzSD2R2p(sd)

% sd: standard deviation of field distribution in units of Gauss, assuming
% a Gaussian distribution

% R2p in units of s-1


gamma=4258*2*pi;

R2p=sd*gamma/sqrt(2);
