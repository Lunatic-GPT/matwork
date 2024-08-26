function y=ff_Gauss_Lorentz(b,x)
% y=ff_Guauss_Lorentz(b,x)
% %4.3f,x_0=%4.3f,f_l=%4.3f,f_g=%4.3f,c=%4.3f

y=b(1).*exp(-(x-b(2))/b(3)).*exp(-(x-b(2)).^2./b(4).^2)+b(5);

