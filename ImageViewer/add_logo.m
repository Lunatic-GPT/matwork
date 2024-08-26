a=imread('LOGO.png');
pvs=imread('pvs.tif');
d=circshift(double(a(:,:,1)'),[0,100]);

pvs2=double(pvs(:,:,1)');
d(268:692,270:end-18)=pvs2;

setappdata(gcf,'udata',d);