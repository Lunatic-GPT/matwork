function R2p=Gz2R2p(Gz,dx,method)

% gz in units of gauss/mm
% dx in units of mm
% R2p in units of s-1
% method: 1 from sinc; 2 from sd and assuming a Gaussian

gamma=4258*2*pi;

if method==1
    x=0.7; %x=fsolve(@(x) sinc(x)-exp(-1),1);
    
    R2p=Gz*gamma*dx/2/x;  % R2p from a sinc function
else
    
    sd=Gz*dx/sqrt(12);
    R2p=BzSD2R2p(sd);
    
end





