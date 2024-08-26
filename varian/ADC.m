function ADC(fname,tind)
% ADC(fname,tind)


a=read_fid(fullfile(fname,'fid'));
%a=read_fid('gdiff_te1p5ms3.fid/fid');
a=squeeze(a);
a=a-repmat(mean(a(end-100:end,:),1),[size(a,1),1]);
gdiff= readPar(fname,'gdiff');

gxd=gdiff;
gyd=gdiff;
gzd=gdiff;

tdefx=readPar(fname,'tdefx');
tdefy=readPar(fname,'tdefy');
tdefz=readPar(fname,'tdefz');
pw90 = readPar(fname,'pw90');
asym = readPar(fname,'asym');

te=readPar(fname,'te');
tm=readPar(fname,'tm');
u2g=40/32767;

tx=tm+2*asym*pw90+tdefx;  % distance between the centers of the diffusion gradients.
ty=tm+te/2-asym*0.26*pw90;
tz=tm+te/2-asym*0.26*pw90;

   bx=(4258*6.28*gxd*u2g*tdefx).*(4258*6.28*gxd*u2g*tdefx)*(tx-tdefx/3)/100;
   by=(4258*6.28*gyd*u2g*tdefy).*(4258*6.28*gyd*u2g*tdefy)*(ty-tdefy/3)/100;
   bz=(4258*6.28*gzd*u2g*tdefz).*(4258*6.28*gzd*u2g*tdefz)*(tz-tdefz/3)/100;
         
bval=bx+by+bz;

aa = abs(a);
b= mean(aa(tind,:),1);

%figure;
hold on;
plot(bval,log10(b),'o');




