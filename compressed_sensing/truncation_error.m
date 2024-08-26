function res=truncation_error(xfm,data,frac)
%truncation_error(xfm,data,frac)
% calculate ||x_0-x_0,s||_l1/sqrt(s)
% second term in Eq. (6) of Candes 2005.
x=xfm*data;

x=sort(abs(x(:)),'Descend');

S=round(frac*length(x));


res=sum(abs(x(S+1:end)))/sqrt(S)/sqrt(length(x));


