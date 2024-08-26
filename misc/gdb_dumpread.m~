%dump memory ref4 ref4.data[0][0] ref4.data[0][0]+128*64*2

sz=[64,64,2];
a=fopen('d','r','ieee-le');

tmp=fread(a,prod(sz)*2,'double');

fclose(a);

re=tmp(1:2:end-1);
im=tmp(2:2:end);

d=re+1i*im;

d=reshape(d,sz);
figure;subplot(1,2,2);imshow(angle(d(:,:,1,1)),[]);
subplot(1,2,1);imshow(abs(d(:,:,1,1)),[]);

tmp2=d(:,:,1,1);
tmp=zrp2(:,:,1,1);

figure;subplot(1,2,2);imshow(angle(tmp(:,:,1,1)),[]);
subplot(1,2,1);imshow(abs(tmp(:,:,1,1)),[]);
