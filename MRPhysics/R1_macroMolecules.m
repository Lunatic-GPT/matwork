% dipolr relaxation constat for macromolecules
%Jokivarsi; jcbfm 2009. 
tau=6.4e-10;
hb=6.626/2/pi*1e-27;
gm=4258*2*pi;
r=1.58e-8;
res=2*1/2*3/2*hb^2*gm^4/r^6*tau; 
disp(res);