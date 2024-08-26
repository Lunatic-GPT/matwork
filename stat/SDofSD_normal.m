function [sdsd,sig]=SDofSD_normal(n,sd)
%Journal of the Royal Statistical Society , 1932, Vol. 95, No. 4 (1932), pp. 695-69
if ~exist('sd','var')
    sd=1;
end

gratio=gamma((n-1)/2)./gamma(n/2);


e_sqrt=(n-1)/2-(1./gratio)^2;


sdsd=sd*gratio.*sqrt(e_sqrt);

sig=sd.*sqrt((n-1)/2).*gratio;