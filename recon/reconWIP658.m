function reconWIP658(fname)
[a,l,s]=readMeasDat(fname);

l=l+1;
s=s+1;
b=reshape(a,[128,2,20,64]);

l=reshape(l,[2,20,64]);
s=reshape(s,[2,20,64]);

c=b;
c(:,:,:,l(1,1,:))=b;
d(:,:,s(1,:,1),:)=c;

e=permute(d,[1,4,3,2]);
nsl=size(e,3);

 if mod(nsl,2)==0
     sorder=[2:2:nsl,1:2:nsl-1];
    else
     sorder=[1:2:nsl-1,2:2:nsl];
    end
    e(:,:,sorder,:)=e;
    
fe=fft2c(e);

prefix=strtok(fname,'.');
writeanalyze(abs(fe),prefix,[1,1,1,1]);

b1=acos(abs(fe(:,:,:,2)./fe(:,:,:,1))-1);

b1=b1/pi*2*500;

writeanalyze(b1,[prefix,'_B1'],[1,1,1,1]);

%figure;imshow(squeeze(abs(fd(:,:,12,1))),[]);

%b=reshape(a,[128,20,2,64]);



