function compare_us_pk(fname,pki,dpi,mask)

ts = BrikLoad(fname);



pk = mean(ts(:,:,:,pki),4);
dp = mean(ts(:,:,:,dpi),4);

x = pk(mask>0);
y = dp(mask>0);
figure;plot(x,y,'ro');

