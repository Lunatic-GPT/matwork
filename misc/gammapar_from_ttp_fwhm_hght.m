function [c,alp,tau,area]=gammapar_from_ttp_fwhm_hght(ttp,fwhm,hght)

%ttp=alp*tau;
tau=fzero(@(x) calc_fwhm(x,ttp/x)-fwhm,ttp*0.8);
alp=ttp/tau;
h1=(tau*alp)^alp*exp(-alp);
c=hght/h1;

area=c*tau^(alp+1)*gamma(alp+1);



function w=calc_fwhm(tau,alp)

x1=fzero(@(x) (tau*alp)^alp*exp(-alp)/2-x^alp*exp(-x/tau),tau*alp*0.5);
x2=fzero(@(x) (tau*alp)^alp*exp(-alp)/2-x^alp*exp(-x/tau),tau*alp*3);
w=x2-x1;